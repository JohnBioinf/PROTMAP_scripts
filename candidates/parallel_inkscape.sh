#!/bin/bash

set -e

data_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['data_dir'])")"

tmp_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['tmp_dir'])")"

selection_cut_off=$1
selection_file="nov_psm$selection_cut_off"

threads=20
pdf_dir="$tmp_dir/candidates_pdf_pics"
pic_dir="$data_dir/candidates/html/$selection_file/pics"
files=$(find $pdf_dir -name "*pdf" -exec basename -s .pdf {} +)
echo "$files" | parallel --tmpdir "$tmp_dir" -j $threads "inkscape -l $pic_dir/{}.svg $pdf_dir/{}.pdf"
