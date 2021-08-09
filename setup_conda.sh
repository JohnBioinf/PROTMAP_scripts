#!/bin/bash

set -e

conda_base=$(conda info | grep -i 'base environment' | awk '{print $4}')
source "$conda_base/etc/profile.d/conda.sh"

conda config --append channels conda-forge
conda config --append channels bioconda
conda config --append channels r

# R
conda create -y -n protmap_R python=3.6

conda activate protmap_R

conda install -y r r-Cairo r-Rcpp r-glue r-ggplot2 r-ggsci r-gridExtra r-htmlwidgets r-latex2exp r-plotly r-RColorBrewer r-rjson r-viridis

conda env export > conda_protmap_R.yml

conda deactivate

# Spectrum Utils
conda create -y -n protmap_spec_util python=3.6

conda activate protmap_spec_util

conda install -y spectrum_utils

conda env export > conda_spectrum_utils.yml

# Rest
conda deactivate

conda create -y -n protmap python=3.6

conda activate protmap

conda install -y blast parallel segemehl ucsc-twobitinfo ucsc-fatotwobit ucsc-bedtobigbed ucsc-wigtobigwig samtools emboss biopython argparse biopython intervaltree jinja2 requests selenium spectrum_utils

pip install splinter

conda env export > conda_protmap.yml
