#!/bin/bash

set -e

echo "################################ Download Genome ################################"
. ./build_db/download_all_genomes.sh

echo "################################ Make Bed ################################"
python ./build_db/gb2bed.py

echo "################################ Build 6frame ################################"
. ./build_db/build_6frame_db.sh

echo "################################ Build Proteom ################################"
python ./build_db/build_proteom_db.py

echo "################################ Download cRAP ################################"
. ./build_db/download_cRAP.sh

echo "################################ Build 6frame To Proteom Map ################################"
python ./build_db/proteom_6frame_map.py

echo "################################ Get validation levels for proteins ################################"
python ./build_db/get_validation_levels.py
