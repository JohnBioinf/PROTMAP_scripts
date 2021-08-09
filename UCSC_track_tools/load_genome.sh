#!/bin/bash

set -e
set -x

genome_dir=$1
track_dir=$2
species=$3

faToTwoBit "$genome_dir/${species}_genome.fasta" "$track_dir/$species.2bit"
twoBitInfo "$track_dir/$species.2bit" stdout | sort -k2rn > "$track_dir/$species.chrom.sizes"
