#!/bin/bash

python gga_codon_muts_oligo_design.py \
    --tiles_csv test_example/KP311_GAA_assembly_fragments.csv \
    --mutations_to_make_csv test_example/mutations_to_make.csv \
    --output_oligos_fasta test_example/output_oligos.fa
