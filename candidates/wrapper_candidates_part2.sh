#!/bin/bash

set -e

selection_cut_off=6
selection_file="nov_psm$selection_cut_off"

echo "Blast candidates"
. ./candidates/parallel_blast.sh "$selection_cut_off"

echo "Collect Results"
python ./candidates/collection.py "$selection_cut_off"

echo "Make HTML"
python ./candidates/make_html.py "$selection_cut_off"

echo "Get data for s-hat"
python ./candidates/scatterplot_s-hat.py "$selection_cut_off"

echo "Build interactiv scatterplots"
Rscript ./candidates/scatterplot_s-hat_html.R "$selection_cut_off"

echo "Get pictures of Candidates"
python ./candidates/make_genome_pics.py "$selection_cut_off"

echo "Convert pictures in parallel"
. ./candidates/parallel_inkscape.sh "$selection_cut_off"
 
echo "Make table for paper"
python ./candidates/make_table_paper.py "$selection_cut_off"
