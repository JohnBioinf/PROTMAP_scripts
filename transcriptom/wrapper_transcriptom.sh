#!/bin/bash

set -e 

segemehl_single() {
	fastq=$1
	genome_dir=$2
	transcriptom_dir=$3

	echo "$fastq"
	segemehl.x -i "$genome_dir/metagenome.idx"\
		-d "$genome_dir/metagenome.fasta"\
		-q "$transcriptom_dir/$fastq.fastq"\
		-u "$transcriptom_dir/$fastq.unpaired_fastq"\
		> "$transcriptom_dir/$fastq.sam"
}

samtools_single() {
	fastq=$1
	genome_dir=$2
	transcriptom_dir=$3

	echo "$fastq"
	samtools view -S -b "$transcriptom_dir/$fastq.sam" > "$transcriptom_dir/$fastq.bam"
	cd "$transcriptom_dir"
	samtools sort -o "${fastq}_sorted.bam" "$transcriptom_dir/$fastq.bam"
	samtools mpileup -S -f "$genome_dir/metagenome.fasta" "${fastq}_sorted.bam" | awk '{print $1, $2-1, $2, $4}' > "$fastq.wig"
}

genome_dir="$data_dir/genome"
transcriptom_dir="$data_dir/transcriptom"
work_dir="$(pwd)"

threads=8

if [ ! -d "$transcriptom_dir" ]; then
	mkdir "$transcriptom_dir"
fi

# get links over https://sra-explorer.info/
ftp_links=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/056/SRR12411056/SRR12411056.fastq.gz"\
	"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/054/SRR12411054/SRR12411054.fastq.gz"\
	"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/057/SRR12411057/SRR12411057.fastq.gz"\
	"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/055/SRR12411055/SRR12411055.fastq.gz"\
	"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/052/SRR12411052/SRR12411052.fastq.gz"\
	"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/059/SRR12411059/SRR12411059.fastq.gz"\
	"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/053/SRR12411053/SRR12411053.fastq.gz"\
	"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR124/058/SRR12411058/SRR12411058.fastq.gz")

i=0
j=0
while [ "$i" -lt $( expr "${#ftp_links[@]}" - 1) ]; do
	if [ "$j" -gt 3 ]; then
		echo "To many tries downloading transcriptom data"
		exit 1
	fi
	ftp_link="${ftp_links[$i]}"
	echo "$ftp_link"
	if wget "$ftp_link" -nv -P "$transcriptom_dir"; then
		file_name=$(basename "$ftp_link")
		gunzip "$transcriptom_dir/$file_name"
		let i=i+1
		j=0
	else
		let j=j+1
	fi
done

fastq_ids=($(find "$transcriptom_dir" -name *.fastq -exec basename -s ".fastq" {} +))

echo "build metagenome"

if [ -f "$genome_dir/metagenome.fasta" ]; then
  rm "$genome_dir/metagenome.fasta"
fi

find "$genome_dir" -name *_genome.fasta -exec cat {} >> "$genome_dir/metagenome.fasta" +

samtools faidx "$genome_dir/metagenome.fasta"

cut -f1-2 "$genome_dir/metagenome.fasta.fai" > "$genome_dir/chrom.sze"
cat "$genome_dir/metagenome.fasta" | grep '>' | awk '{id = substr($1, 2); name = $2 " " $3; print id ";" name}' > "$genome_dir/metagenome_header.hash" 

echo "segemehl"

# build segemehl index
segemehl.x --silent -x "$genome_dir/metagenome.idx" -d "$genome_dir/metagenome.fasta" -t $threads

export -f segemehl_single
parallel --tmpdir "$tmp_dir" -j "$threads" segemehl_single {} "$genome_dir" "$transcriptom_dir" ::: "${fastq_ids[@]}"

echo "coverage"

export -f samtools_single
parallel --tmpdir "$tmp_dir" -j "$threads" samtools_single {} "$genome_dir" "$transcriptom_dir" ::: "${fastq_ids[@]}"

echo "build average"

python ./transcriptom/make_average.py 

cd "$transcriptom_dir"
wigToBigWig average.wig "$genome_dir/chrom.sze" average.bw
cd "$work_dir"
