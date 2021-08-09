#!/bin/bash

set -e

cand_dir="$data_dir/candidates"

selection_cut_off=6
selection_file="nov_psm$selection_cut_off"

bed_dir="$cand_dir/bed_dir/$selection_file"
blast_dir="$cand_dir/blast/$selection_file"
html_dir="$cand_dir/html/$selection_file"

if [ ! -d "$cand_dir" ]; then
	mkdir "$cand_dir"
fi

if [ ! -d "$bed_dir" ]; then
	mkdir -p "$bed_dir"
fi

if [ ! -d "$blast_dir" ]; then
	mkdir -p "$blast_dir"
fi

if [ ! -d "$html_dir" ]; then
	mkdir -p "$html_dir/pics"
fi

echo "Select candidates"
python ./candidates/selection.py "$selection_cut_off"

echo "Make bed files for candidates"
python ./candidates/make_bed_and_fasta.py "$selection_cut_off"
