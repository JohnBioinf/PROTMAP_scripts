#!/bin/bash

set -e
set -x
species=$1
output_dir=$2
massSpec_dir=$3
db=$4
tmp_dir=$5
organism=$6

master_massSpec_dir="$tmp_dir/$db/$species"

if [ -d "$master_massSpec_dir" ]; then
	rm "$master_massSpec_dir/"*
else
	mkdir -p "$master_massSpec_dir"
fi

if [ -f "$master_massSpec_dir/all_unsorted.bed" ]; then
	rm "$master_massSpec_dir/all_unsorted.bed"
fi

if [ -f "$master_massSpec_dir/all.bed" ]; then
	rm "$master_massSpec_dir/all.bed"
fi

sed '/track name=.*/d' "$massSpec_dir/$species/"*".bed" >> "$master_massSpec_dir/all_unsorted.bed"
printf "\n" >> "$master_massSpec_dir/all_unsorted.bed"

sort -k1,1 -k2,2n "$master_massSpec_dir/all_unsorted.bed" > "$master_massSpec_dir/all.bed"
bedToBigBed -extraIndex=name "$master_massSpec_dir/all.bed" "$output_dir/$species.chrom.sizes" "$output_dir/bbi/massSpec_${db}_${organism}.bb"
