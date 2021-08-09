#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

input_file=$1
current_db=$2
res_dir=$3
protein_db_type=$4

project_name=$(echo $input_file | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
echo $project_name
if [[ $project_name == *"AspN"* ]]; then
	param_file="$(pwd)/comet/comet_${protein_db_type}_AspN.params"
else
	param_file="$(pwd)/comet/comet_$protein_db_type.params"
fi

comet.2019014.linux.exe "$input_file" -N"$res_dir/$project_name" -D"$current_db" -P$param_file
#awk 'BEGIN{i = 0}{print i++ "\t" $0}' "$res_dir/${project_name}.txt" > "$res_dir/${project_name}_id.txt"
#awk 'BEGIN{i = 0}{print i++ "\t" $0}' "$res_dir/${project_name}.decoy.txt" > "$res_dir/${project_name}_id.decoy.txt"
