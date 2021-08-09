#!/bin/bash

set -e

db_dir="$data_dir/dbs"

if [ -f "$db_dir/crap.fasta" ]; then
	rm "$db_dir/crap.fasta"
fi

cRAP_link="ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta"

wget "$cRAP_link" -P "$db_dir"

cat "$db_dir/SIHUMI_6frame.fasta" "$db_dir/crap.fasta" > "$db_dir/SIHUMI_6frame_cRAP.fasta"
cat "$db_dir/ecoli_6frame.fasta" "$db_dir/crap.fasta" > "$db_dir/ecoli_6frame_cRAP.fasta"

cat "$db_dir/SIHUMI_proteom.fasta" "$db_dir/crap.fasta" > "$db_dir/SIHUMI_proteom_cRAP.fasta"
cat "$db_dir/ecoli_proteom.fasta" "$db_dir/crap.fasta" > "$db_dir/ecoli_proteom_cRAP.fasta"
