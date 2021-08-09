#!/bin/bash
# set -x
set -e

species=$1
output_dir=$2
tmp_dir=$3
genome_dir=$4

conda_base=$(conda info | grep -i 'base environment' | awk '{print $4}')
source "$conda_base/etc/profile.d/conda.sh"
conda activate protmap

python ./UCSC_track_tools/gb2bed.py "$species" "$genome_dir" "/$tmp_dir/$species.bed"
sed '/track name=.*/d' "/$tmp_dir/$species.bed" > "/$tmp_dir/${species}_unsorted.bed"
sort -k1,1 -k2,2n "/$tmp_dir/${species}_unsorted.bed" > "/$tmp_dir/$species.bed"
bedToBigBed -extraIndex=name "/$tmp_dir/$species.bed" "$output_dir/$species.chrom.sizes" "$output_dir/bbi/annotation.bb"
