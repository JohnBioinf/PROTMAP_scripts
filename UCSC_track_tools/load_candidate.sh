#!/bin/bash

set -x
set -e

cand_dir=$1
species=$2
output_dir=$3
tmp_dir=$4

sed '/track name=.*/d' "/$cand_dir/${species}_cand.bed" > "/$tmp_dir/${species}_unsorted.bed"
sort -k1,1 -k2,2n "/$tmp_dir/${species}_unsorted.bed" > "/$tmp_dir/$species.bed"
bedToBigBed -extraIndex=name "/$tmp_dir/$species.bed" "$output_dir/$species.chrom.sizes" "$output_dir/bbi/cand.bb"
