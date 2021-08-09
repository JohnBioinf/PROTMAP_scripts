#!/bin/bash

set -e

selection_cut_off=$1
selection_file="nov_psm$selection_cut_off"
echo "$selection_file"

data_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['data_dir'])")"
tmp_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['tmp_dir'])")"

blastdb_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['blastdb_dir'])")"
blast_dir="$data_dir/candidates/blast/$selection_file"

threads=20
ls "$blast_dir"/*fasta | parallel --tmpdir "$tmp_dir" -j $threads --joblog "./candidates/blast_parallel.log" "blastp -db $blastdb_dir -query {} -outfmt '6 staxids evalue sseq pident sacc salltitles ssciname qcovs' -out {}.tsv"
