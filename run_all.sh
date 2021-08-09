#!/bin/bash

set -e

conda_base=$(conda info | grep -i 'base environment' | awk '{print $4}')
source "$conda_base/etc/profile.d/conda.sh"
conda activate protmap

data_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['data_dir'])")"
tmp_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['tmp_dir'])")"
bin_path="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['bin_path'])")"
ms_dir="$(cat "./parameters.json" |
	python3 -c "import sys, json; print(json.load(sys.stdin)['ms_dir'])")"

PATH="$bin_path:$PATH"

# 1. Download genomes and build data bases
echo "Download and build all necesary data bases"
. ./build_db/wrapper_build_db.sh

# 2. Download transcriptom data and map reads
echo "Process transcriptomic data"
. ./transcriptom/wrapper_transcriptom.sh

# 3. Download MS-Spectra and build peptide to spectra map
echo "comet"
. ./comet/wrapper_comet.sh

# 4. Aggregate PSMs
echo "Aggregate PSMs"
python ./data_accumulation/aggreagate_psm.py

# 5. Candidates Part 1 (Select all candidates and make bed)
echo "Select candidates"
. ./candidates/wrapper_candidates_part1.sh

# 6. Build UCSC Track Hub
echo "Build Trackhub"
python ./UCSC_track_tools/build_mainstructure.py

# Might put a break here to upload UCSCS Track Hub if track hub is not yet loaded
echo "Upload Track hub to UCSC. Edit parameters.json accordingly. Resume"
# exit 0

# 7. Candidates Part 2 (Make html result page)
echo "Build candidates result page"
. ./candidates/wrapper_candidates_part2.sh

# 8. Generate all data for plotting
echo "Collect data for Plots"
. ./data_accumulation/wrapper_accumulation.sh

# 9. Figure Plotting
echo "Plot figures"
. ./figure_plotting/wrapper_plotting.sh

# 10. Wrong start annotation ecoli html list.
echo "Build html for ecoli start annotation errors."
python ./start_anno_html/wrapper_start_anno_html.py

conda activate spec_util
python ./start_anno_html/wrapper_start_anno_html_2.py

echo "Protmap pipeline finished"
