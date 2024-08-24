# Script for designing oligos for Golden-Gate Assembly codon mutagenesis

This Python script can be used to design oligos for Golden-Gate Assembly
codon mutagenesis. It requires you to create a CSV containing the tiles
fragments and a CSV specifying the mutations to make.
Those are then passed to the script as arguments (see the command-line usage below),
and the script can be run to generate the oligos. You can specify different
representations for mutations (if you want some mutations covered by multiple oligos),
as well as various other options including motifs (eg, restriction sites) to avoid;
this is all done using the command-line arguments shown in the usage below.

To run it, just download [gga_codon_muts_oligo_design.py](gga_codon_muts_oligo_design.py)
and then run with:
```
python gga_codon_muts_oligo_design.py <arguments>
```

If you want to run the [test example](test_example), clone this entire repository and
then run the command in [run_on_test_example.bash](run_on_test_example.bash), eg:
```
python gga_codon_muts_oligo_design.py \
    --tiles_csv test_example/KP311_GAA_assembly_fragments.csv \
    --mutations_to_make_csv test_example/mutations_to_make.csv \
    --output_oligos_fasta test_example/output_oligos.fa
```

Note that there other tools that do similar (and for some uses, possibly more or better)
versions of the same tasks implemented in this script, such as
[DIMPLE](https://github.com/coywil26/DIMPLE).

## Requirements
This script requires that you have installed:
 - Python >= version 3.8
 - A relatively recent version of BioPython
 - A relatively recent version of pandas


## Usage

The command-line usage of the script can be seen by running the command
`python gga_codon_muts_oligo_design.py`, and is shown below:

```
usage: gga_codon_muts_oligo_design.py [-h] --tiles_csv TILES_CSV --mutations_to_make_csv
                                      MUTATIONS_TO_MAKE_CSV --output_oligos_fasta
                                      OUTPUT_OLIGOS_FASTA
                                      [--max_representation MAX_REPRESENTATION]
                                      [--wildtype_frac WILDTYPE_FRAC]
                                      [--avoid_motifs AVOID_MOTIFS]
                                      [--codon_freqs_csv CODON_FREQS_CSV]

Design oligos for tiles for Golden-Gate assembly codon mutagenesis. To use this script,
first you need to break your gene into tiles of that can be ordered (be sure to design tiles
that will give good overhangs; https://pubs.acs.org/doi/10.1021/acssynbio.8b00333). You then
specify those tiles using the '--tiles_csv' argument, and also specify the mutations to make
and the representation (number of oligos) for each one in '--mutations_to_make_csv'. A
representation of 1 means a single oligo for that mutation is made; larger representation
values mean more oligos for each mutation are made which should increase its representation
in the final library. See also '--max_representation'. If multiple oligos are made for the
same mutation, when possible they use different codons.

options:
  -h, --help            show this help message and exit
  --tiles_csv TILES_CSV
                        CSV with nucleotide sequences of tiles, should have columns
                        'fragment' (fragment name), 'fragment_sequence' (full nucleotide
                        sequence of fragment), and 'inframe_mutated_region' (nucleotide
                        sequence of part of fragment that is in-frame mutated region of
                        gene). Fragments must be in order that their
                        'inframe_mutation_region' sequences should be concatenated to make
                        the full gene. The overall 'fragment_sequence' will have flanking
                        regions for Golden Gate assembly that are not present in
                        'inframe_mutated_region'. Be sure to specify any restriction enzymes
                        that will be used in '--avoid_motifs'. (default: None)
  --mutations_to_make_csv MUTATIONS_TO_MAKE_CSV
                        CSV with mutations to make. Must include columns 'sequential_site'
                        (site number in 1, 2, numbering of protein), 'wildtype_aa' (parental
                        amino acid at that site), 'mutant_aa' (the mutation to make at the
                        site), and 'representation' (how many oligos to make with that
                        mutation; see also '--max-representation'). (default: None)
  --output_oligos_fasta OUTPUT_OLIGOS_FASTA
                        Output FASTA file with created oligos. The oligos are named
                        according to the sequential site that is mutated (not the reference
                        site) (default: None)
  --max_representation MAX_REPRESENTATION
                        The maximum representation (number of oligos) for any mutation
                        regardless of value given in '--mutations_to_make_csv'. (default: 2)
  --wildtype_frac WILDTYPE_FRAC
                        For each tile, a wildtype sequence is included to an amount equal to
                        ceiling of this fraction times the number of mutations for that
                        tile. (default: 0.005)
  --avoid_motifs AVOID_MOTIFS
                        Avoid these motifs and reverse complements (typically restrition
                        sites). (default: ['CGTCTC'])
  --codon_freqs_csv CODON_FREQS_CSV
                        File specifying a frequency for each codon for an amino acid. Codons
                        are chose to first prioritize the highest-frequency one for that
                        amino acid. Must have columns 'codon', 'aa', and 'frequency'.
                        (default: https://raw.githubusercontent.com/jbloomlab/gga_codon_muts
                        _oligo_design/main/human_codon_freq.csv)
```
