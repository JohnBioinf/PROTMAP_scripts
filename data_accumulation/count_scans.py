#!/usr/bin/python3

import os
import subprocess
from pyteomics import mzml

mzml_files = subprocess.check_output("find /scr/k61san2/john/PROTMAP_data/ms/ -name *.mzML", shell=True).decode().split("\n")

sum_scans = 0
for mzml_file in mzml_files:
    if not os.path.isfile(mzml_file):
        continue
    with mzml.read(mzml_file) as reader:
        for scan in reader:
            sum_scans += 1


with open("../spectral_count.txt", "w") as file_handle:
    file_handle.write("The number of spectras over all experiments are: "
                      + str(sum_scans) + ".\n")
