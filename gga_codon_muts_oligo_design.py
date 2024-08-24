#!/usr/bin/env python3

"""Script to design oligos for Golden-Gate assembly codon mutagenesis."""


import argparse
import itertools
import math
import re
import sys

import Bio.SeqIO

import pandas as pd


MIN_PYTHON_VERSION = (3, 8)
if sys.version_info < MIN_PYTHON_VERSION:
    raise RuntimeError(
        f"Script requires Python >= {MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]}"
    )


def remove_motif(oligo, motif, aa_to_codon):
    """Remove motif from oligo while keeping protein sequence."""
    assert len(oligo) % 3 == 0
    if motif not in oligo:
        return oligo
    prot = str(Bio.Seq.Seq(oligo).translate())

    while motif in oligo:
        i = oligo.index(motif)
        for icodon in range(i // 3, i // 3 + 1 + len(motif) // 3):
            aa = prot[icodon]
            codon = oligo[icodon * 3: icodon * 3 + 3]
            other_codons = [c for c in aa_to_codon[aa] if c != codon]
            for other_codon in other_codons:
                oligo = oligo[: icodon * 3] + other_codon + oligo[icodon * 3 + 3: ]
                if oligo[i: i + len(motif)] != motif:
                    break

            if oligo[i: i + len(motif)] != motif:
                break

        if oligo[i: i + len(motif)] == motif:
            raise ValueError(f"Cannot remove {motif=} from {oligo=} at {i=}")

        if (motif in oligo) and (oligo.index(motif) <= i):
            raise ValueError(
                f"Removing {motif=} at {i=} from {oligo=} created earlier motif"
            )

    assert motif not in oligo
    assert prot == str(Bio.Seq.Seq(oligo).translate())

    return oligo


def gga_codon_muts_oligo_design(
    tiles_csv,
    mutations_to_make_csv,
    output_oligos_fasta,
    max_representation,
    wildtype_frac,
    avoid_motifs,
    codon_freqs_csv,
):
    """Function that implements the oligo design."""
    print(f"\nReading tiles from {tiles_csv=}")
    tiles = (
        pd.read_csv(tiles_csv)
        [["fragment", "fragment_sequence", "inframe_mutated_region"]]
        .assign(
            fragment_sequence=lambda x: x["fragment_sequence"].str.upper(),
            inframe_mutated_region=lambda x: x["inframe_mutated_region"].str.upper(),
            inframe_start=lambda x: x.apply(
                lambda row: row["fragment_sequence"].find(row["inframe_mutated_region"]),
                axis=1,
            ),
            upstream_flank=lambda x: x.apply(
                lambda row: row["fragment_sequence"][0: row["inframe_start"]],
                axis=1,
            ),
            downstream_flank=lambda x: x.apply(
                lambda row: row["fragment_sequence"][
                    row["inframe_start"] + len(row["inframe_mutated_region"]):
                ],
                axis=1,
            ),
        )
    )
    assert len(tiles) == tiles["fragment"].nunique(), "tiles not uniquely named"

    records = []
    prot_from_tiles = []
    sequential_start = 1
    for tup in tiles.itertuples():
        fragment_prot = str(Bio.Seq.Seq(tup.inframe_mutated_region).translate())
        if len(tup.inframe_mutated_region) % 3 != 0:
            raise ValueError(
                f"'inframe_mutated_region' of fragment {tup.fragment} has length "
                "that is not a multiple of 3:\n"
                f"length: {len(tup.inframe_mutated_region)}\n"
                f"inframe_mutated_region: {tup.inframe_mutated_region}"
            )
        if "*" in tup:
            raise ValueError(
                f"in-frame region of fragment {tup.fragment} encodes a stop codon"
            )
        if "-" in tup:
            raise ValueError(
                f"in-frame region of fragment {tup.fragment} encodes a gap"
            )
        prot_from_tiles.append(fragment_prot)
        sequential_end = sequential_start + len(fragment_prot) - 1
        records.append((tup.fragment, fragment_prot, sequential_start, sequential_end))
        sequential_start = sequential_end + 1
    prot_from_tiles = "".join(prot_from_tiles)
    assert len(prot_from_tiles) == sequential_end, f"{len(prot_from_tiles)=}, {sequential_end=}"
    print(f"Tiles encode protein of {len(prot_from_tiles)} residues:\n{prot_from_tiles}\n")
    tiles = tiles.merge(
        pd.DataFrame(
            records,
            columns=["fragment", "fragment_prot", "sequential_start", "sequential_end"]
        ),
        on="fragment",
        validate="one_to_one",
    )

    print(f"Reading mutations to make from {mutations_to_make_csv=}")
    mutations_to_make = pd.read_csv(mutations_to_make_csv)[
        ["sequential_site", "wildtype_aa", "mutant_aa", "representation"]
    ]
    assert len(mutations_to_make) == len(mutations_to_make.drop_duplicates())
    if any(mutations_to_make["wildtype_aa"] == mutations_to_make["mutant_aa"]):
        raise ValueError(
            "'mutations_to_make_csv' has some sites where 'wildtype_aa' = 'mutant_aa'; "
            "remove these."
        )
    if len(mutations_to_make) != len(
        mutations_to_make.groupby(["sequential_site", "mutant_aa"])
    ):
        raise ValueError(
            "Rows in 'mutations_to_make_csv' do not each specify unique "
            "'sequential_site' and 'mutant_aa'."
        )

    # check protein in `mutations_to_make` matches that encoded by tiles
    prot_to_make = mutations_to_make[["sequential_site", "wildtype_aa"]].drop_duplicates()
    if len(prot_to_make) != prot_to_make["sequential_site"].nunique():
        raise ValueError(
            "'mutations_to_make_csv' has multiple 'wildtype_aa' for some 'sequential_site'"
        )
    prot_to_make = prot_to_make.set_index("sequential_site")["wildtype_aa"].to_dict()
    for r, aa in prot_to_make.items():
        if r > len(prot_from_tiles):
            raise ValueError(
                f"'sequential_site' {r} in 'mutations_to_make_csv' it outside range of "
                "protein specified in 'tiles_csv'"
            )
        if prot_from_tiles[r - 1] != aa:
            raise ValueError(
                f"At 'sequential_site' {r}, mismatch in 'wildtype_aa' in "
                "'mutations_to_make_csv' and protein encoded in 'tiles_csv': "
                f"{aa} versus {prot_from_tiles[r - 1]}"
            )
    print("Representation values for the mutations to make:")
    print(
        mutations_to_make
        .groupby("representation")
        .aggregate(n_mutations=pd.NamedAgg("sequential_site", "count"))
    )
    mutations_to_make = mutations_to_make.query("representation > 0").assign(
        representation=lambda x: x["representation"].clip(upper=max_representation)
    )
    print("Representations after removing zeros and clipping at {max_representation=}")
    print(
        mutations_to_make
        .groupby("representation")
        .aggregate(n_mutations=pd.NamedAgg("sequential_site", "count"))
    )
    print(
        f"So overall, we will make {len(mutations_to_make)} mutations encompassing "
        f"{mutations_to_make['representation'].sum()} non-wildtype oligos.\n"
    )

    avoid_motifs = set(
        [m.upper() for m in avoid_motifs]
        + [str(Bio.Seq.Seq(m).reverse_complement()) for m in avoid_motifs]
    )
    print(f"We will avoid the following motifs:\n{avoid_motifs}\n")
    assert all(re.fullmatch("[ATCG]+", m) for m in avoid_motifs)
    nt_from_tiles = "".join(tiles["inframe_mutated_region"])
    for motif in avoid_motifs:
        if motif in nt_from_tiles:
            raise ValueError(
                f"{motif=} is already in parent nucleotide sequence in 'tiles_csv'"
            )

    print(f"Reading the codon frequencies to use from {codon_freqs_csv=}\n")
    codon_freqs = pd.read_csv(codon_freqs_csv)[["codon", "aa", "frequency"]].assign(
        codon=lambda x: x["codon"].str.upper(),
        aa=lambda x: x["aa"].str.upper(),
    )
    assert len(codon_freqs) == 64
    possible_codons = {
        "".join(tup) for tup in itertools.product(["A", "C", "T", "G"], repeat=3)
    }
    assert set(codon_freqs["codon"]) == possible_codons, f"{codon_freqs['codon']=}\n{possible_codons=}"
    aa_to_codon = (
        codon_freqs
        .sort_values("frequency", ascending=False)
        .groupby("aa")
        .aggregate(codons=pd.NamedAgg("codon", list))
        ["codons"]
        .to_dict()
    )
    aa_to_codon["-"] = [""]
    assert len(aa_to_codon) == 22, f"{len(aa_to_codon)=}"

    # design the oligos
    oligos = []
    n_mut_oligos = 0
    for tile_tup in tiles.itertuples():
        fragment = tile_tup.fragment
        start = tile_tup.sequential_start
        end = tile_tup.sequential_end
        upstream_flank = tile_tup.upstream_flank.lower()
        downstream_flank = tile_tup.downstream_flank.lower()
        ntseq_by_codon = [
            tile_tup.inframe_mutated_region[3 * r: 3 * r + 3]
            for r in range(len(tile_tup.inframe_mutated_region) // 3)
        ]
        assert "".join(ntseq_by_codon) == tile_tup.inframe_mutated_region
        tile_muts = mutations_to_make.query(
            "(sequential_site >= @start) and (sequential_site <= @end)"
        )
        n_tile_mut_oligos = tile_muts["representation"].sum()
        n_tile_wt_oligos = int(math.ceil(n_tile_mut_oligos * wildtype_frac))
        print(
            f"For tile {fragment=} making {len(tile_muts)} mutations with "
            f"{n_tile_mut_oligos} oligos; also {n_tile_wt_oligos} wildtype oligos."
        )
        if len(tile_muts) == 0:
            raise ValueError(f"No mutations to make for tile {fragment=}")
        oligos += [
            (
                f"tile-{fragment}_wildtype_{i + 1}",
                upstream_flank + "".join(ntseq_by_codon) + downstream_flank
            )
            for i in range(n_tile_wt_oligos)
        ]
        for mut_tup in tile_muts.itertuples():
            r = mut_tup.sequential_site - start
            wt_codon = ntseq_by_codon[r]
            wt_aa = str(Bio.Seq.Seq(wt_codon).translate())
            assert mut_tup.wildtype_aa == wt_aa, f"{wt_aa=}, {wt_codon=}, {mut_tup.wildtype_aa=}"
            for i, mut_codon in zip(
                range(mut_tup.representation), itertools.cycle(aa_to_codon[mut_tup.mutant_aa])
            ):
                oligo_name = f"tile-{fragment}_{wt_aa}{mut_tup.sequential_site}{mut_tup.mutant_aa}_{i + 1}"
                oligo = "".join(ntseq_by_codon[: r] + [mut_codon] + ntseq_by_codon[r + 1:])
                for motif in avoid_motifs:
                    if motif in oligo:
                        oligo = remove_motif(oligo, motif, aa_to_codon)
                oligos.append((oligo_name, upstream_flank + oligo + downstream_flank))
                n_mut_oligos += 1
    assert n_mut_oligos == mutations_to_make["representation"].sum() <= len(oligos)
    print(f"\nOverall designed {len(oligos)} oligos including the wildtype ones.")
    nunique = len(set(tup[1] for tup in oligos))
    print(f"{nunique} of these oligos have unique sequences.")

    print(f"\nWriting the oligos to {output_oligos_fasta=}")
    with open(output_oligos_fasta, "w") as f:
        f.write("".join(f">{oligo_name}\n{oligo}\n" for (oligo_name, oligo) in oligos))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Design oligos for tiles for Golden-Gate assembly codon mutagenesis. "
            "To use this script, first you need to break your gene into tiles of "
            "that can be ordered (be sure to design tiles that will give good "
            "overhangs; https://pubs.acs.org/doi/10.1021/acssynbio.8b00333). "
            "You then specify those tiles using the '--tiles_csv' argument, and "
            "also specify the mutations to make and the representation (number of "
            "oligos) for each one in '--mutations_to_make_csv'. A representation of 1 "
            "means a single oligo for that mutation is made; larger representation "
            "values mean more oligos for each mutation are made which should increase "
            "its representation in the final library. See also '--max_representation'. "
            "If multiple oligos are made for the same mutation, when possible they "
            "use different codons."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--tiles_csv",
        help=(
            "CSV with nucleotide sequences of tiles, should have columns 'fragment' "
            "(fragment name), 'fragment_sequence' (full nucleotide sequence of fragment)"
            ", and 'inframe_mutated_region' (nucleotide sequence of part of fragment "
            "that is in-frame mutated region of gene). Fragments must be in order that "
            "their 'inframe_mutation_region' sequences should be concatenated to make "
            "the full gene. The overall 'fragment_sequence' will have flanking regions "
            " for Golden Gate assembly that are not present in 'inframe_mutated_region'. "
            "Be sure to specify any restriction enzymes that will be "
            "used in '--avoid_motifs'."
        ),
        required=True,
    )

    parser.add_argument(
        "--mutations_to_make_csv",
        help=(
            "CSV with mutations to make. Must include columns 'sequential_site' ("
            "site number in 1, 2, numbering of protein), 'wildtype_aa' (parental "
            "amino acid at that site), 'mutant_aa' (the mutation to make at the "
            "site), and 'representation' (how many oligos to make with that "
            "mutation; see also '--max-representation')."
        ),
        required=True,
    )

    parser.add_argument(
        "--output_oligos_fasta",
        help=(
            "Output FASTA file with created oligos. The oligos are named according "
            "to the sequential site that is mutated (not the reference site)"
        ),
        required=True,
    )

    parser.add_argument(
        "--max_representation",
        default=2,
        help=(
            "The maximum representation (number of oligos) for any mutation "
            "regardless of value given in '--mutations_to_make_csv'."
        ),
        type=int,
    )

    parser.add_argument(
        "--wildtype_frac",
        default=0.005,
        help=(
            "For each tile, a wildtype sequence is included to an amount equal to "
            "ceiling of this fraction times the number of mutations for that tile."
        ),
        type=float,
    )

    parser.add_argument(
        "--avoid_motifs",
        help="Avoid these motifs and reverse complements (typically restrition sites).",
        default=["CGTCTC"],
        nargs="+",
    )

    parser.add_argument(
        "--codon_freqs_csv",
        help=(
            "File specifying a frequency for each codon for an amino acid. Codons are "
            "chose to first prioritize the highest-frequency one for that amino acid. "
            "Must have columns 'codon', 'aa', and 'frequency'."
        ),
        default="https://raw.githubusercontent.com/jbloomlab/gga_codon_muts_oligo_design/main/human_codon_freq.csv",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    args = parser.parse_args()

    gga_codon_muts_oligo_design(**vars(args))
