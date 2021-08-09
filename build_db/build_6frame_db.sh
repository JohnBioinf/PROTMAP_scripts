#!/bin/bash

data_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['data_dir'])")"
db_dir="$data_dir/dbs"
if ! [ -d "$db_dir" ]; then
       mkdir "$db_dir"
fi

SIHUMI_db_file="$data_dir/dbs/SIHUMI_6frame.fasta"
ecoli_db_file="$data_dir/dbs/ecoli_6frame.fasta"

genome_dir="$data_dir/genome"
if [ -f "$SIHUMI_db_file" ]; then
	rm $SIHUMI_db_file
fi

if [ -f "$ecoli_db_file" ]; then
	rm $ecoli_db_file
fi

for file in "$genome_dir/"*.gbk; do
	spec=$(basename $file .gbk)
	genome_file="$genome_dir/$(basename $file .gbk)_genome.fasta"
	frame_file="$genome_dir/$(basename $file .gbk)_6frame.fasta"

	if ! [ -f $genome_file ]; then
		python ./build_db/gbk_to_fasta.py "$file" "$genome_file"
	fi
	getorf -table 11 -sequence $file -outseq $frame_file.temp -reverse -minsize 3
	contigs=($(awk '{if($0 ~ ">"){print substr($1,2)}}' "$genome_file"))
	awk -v spec="$spec" -v con="${contigs[*]}" 'BEGIN{
	split(con,contigs," ")
	}{
	if($0 ~ ">"){\
	       	if($0 ~ "REVERSE SENSE"){\
			frame = "-1";\
			start = substr($4,1,length($4) - 1) - 1;\
			stop = substr($2,2)}\
		else{\
			frame = "1";\
			stop = substr($4,1,length($4) - 1);\
			start = substr($2,2) - 1}\
		n = split($1, ident_array, "_");\
		contig = substr(ident_array[1], 2)
		for (i=1;i<=length(contigs);i++){
			if(contigs[i] ~ contig){
				contig = contigs[i]
			}
		}
		print ">" spec "|" contig "|" ident_array[n] "|" frame "|" start "-" stop}\
	else{
		print $0}}' $frame_file.temp > $frame_file
	rm $frame_file.temp
	cat $frame_file >> $SIHUMI_db_file
	if [[ "$spec" == "ecoli" ]]; then
		cat $frame_file >> $ecoli_db_file
	fi
done
