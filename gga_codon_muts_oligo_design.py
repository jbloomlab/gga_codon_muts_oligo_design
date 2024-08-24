#!/usr/bin/env python3

"""Script to design oligos for Golden-Gate assembly codon mutagenesis."""


import argparse
import re
import sys

import Bio.SeqIO

import pandas as pd


MIN_PYTHON_VERSION = (3, 8)
if sys.version_info < MIN_PYTHON_VERSION:
    raise RuntimeError(
        f"Script requires Python >= {MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]}"
    )


def gga_codon_muts_oligo_design(
    tiles_csv,
    mutations_to_make_csv,
    output_oligos_fasta,
    max_representation,
    wildtype_frac,
    avoid_motifs,
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
            upstream_flank=lambda x: x["fragment_sequence"].str.slice(
                start=0, stop=x["inframe_start"]
            ),
            downstream_flank=lambda x: x["fragment_sequence"].str.slice(
                start=x["inframe_start"] + x["inframe_mutated_region"].map(len)
            )
        )
        .drop(columns=["inframe_start"])
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
    print(f"We will avoid the following motifs:\n{avoid_motifs}")
    assert all(re.fullmatch("[ATCG]+", m) for m in avoid_motifs)
    nt_from_tiles = "".join(tiles["inframe_mutated_region"])
    for motif in avoid_motifs:
        if motif in nt_from_tiles:
            raise ValueError(
                f"{motif=} is already in parent nucleotide sequence in 'tiles_csv'"
            )

    raise NotImplementedError
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Design oligos for tiles for Golden-Gate assembly codon mutagenesis. "
            "To use this script, first you need to break your gene into tiles of "
            "that can be ordered. (be sure to design tiles that will give good "
            "overhangs; https://pubs.acs.org/doi/10.1021/acssynbio.8b00333). "
            "You then specify those tiles using the '--tiles_csv' argument, and "
            "also specify the mutations to make and the representation (number of "
            "oligos) for each one in '--mutations_to_make_csv'. A representation of 1 "
            "means a single oligo for that mutation is made; larger representation "
            "values mean more oligos for each mutation are made which should increase "
            "its representation in the final library. See also '--max_representation'."
        )
    )

    parser.add_argument(
        "--tiles_csv",
        help=(
            "CSV with nucleotide sequences of tiles, should have columns 'fragment' "
            "(fragment name), 'fragment_sequence' (full nucleotide sequence of fragment)"
            ", and 'inframe_mutated_region' (nucleotide sequence of part of fragment ",
            "that is in-frame mutated region of gene. Fragments must be in order that "
            "their 'inframe_mutation_region' sequences should be concatenated to make "
            "the full gene. Be sure to specify any restriction enzymes that will be "
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
        help="Output FASTA file with created oligos.",
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
        default=0.001,
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
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    args = parser.parse_args()

    gga_codon_muts_oligo_design(**vars(args))
