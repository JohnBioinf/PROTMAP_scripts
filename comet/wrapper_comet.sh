#!/bin/bash

set -e

db_dir="$data_dir/dbs"
ms_dir="$data_dir/ms"
threads=20

if [ ! -d "$ms_dir" ]; then
	mkdir -p "$ms_dir/SIHUMI/mzML"
	mkdir -p "$ms_dir/ecoli/mzML"

	mkdir -p "$ms_dir/SIHUMI/psm/6frame"
	mkdir -p "$ms_dir/ecoli/psm/6frame"
	mkdir "$ms_dir/SIHUMI/psm/proteom"
	mkdir "$ms_dir/ecoli/psm/proteom"

	mkdir -p "$ms_dir/SIHUMI/bed/6frame"
	mkdir -p "$ms_dir/ecoli/bed/6frame"
	mkdir "$ms_dir/SIHUMI/bed/proteom"
	mkdir "$ms_dir/ecoli/bed/proteom"
fi

# This is temporarily until paper got published and PRIDE data is publicly available

find "$ms_dir" -name "*ecoli*mzML" -exec ln -s {} "$ms_dir/ecoli/mzML/" \;

find "$ms_dir" -name "*SIHUMI*mzML" -exec ln -s {} "$ms_dir/SIHUMI/mzML/" \;

mzML_files=($(find $ms_dir -name "*.mzML"))
db_types=("6frame" "proteom")

touch ./comet/parallel_comet.commands
rm ./comet/parallel_comet.commands

for db_type in ${db_types[@]}; do
	for mzML_file in ${mzML_files[@]}; do
		psm_dir="${mzML_file/\/mzML\//\/psm\/}"
		species="$(echo "$psm_dir" | awk -F '/' '{print $(NF-2)}')"
		psm_dir="$(dirname "${psm_dir}")/$db_type"
		current_db="$db_dir/${species}_${db_type}_cRAP.fasta"
		echo $mzML_file $current_db $psm_dir $db_type >> ./comet/parallel_comet.commands
	done
done

# calling comet in parallel
cat ./comet/parallel_comet.commands | parallel --tmpdir "$tmp_dir" --colsep ' ' -j $threads --joblog "./comet/comet_parallel.log" "./comet/single_comet.sh {1} {2} {3} {4}"

touch ./comet/parallel_pepG.commands
rm ./comet/parallel_pepG.commands

psm_files=($(find "$ms_dir" -name "*.txt"))
for psm_file in ${psm_files[@]}; do
	db_type="$(echo "$psm_file" | awk -F '/' '{print $(NF-1)}')"
	species="$(echo "$psm_file" | awk -F '/' '{print $(NF-3)}')"
	res_dir="$ms_dir/$species/bed/$db_type"
	echo "$psm_file $res_dir $db_type" >> ./comet/parallel_pepG.commands
done

cat ./comet/parallel_pepG.commands | parallel --tmpdir "$tmp_dir" --colsep ' ' -j "$threads" --joblog "./comet/pepG_parallel.log" "python ./comet/pepG_v1.2.py {1} {2} {3}"
