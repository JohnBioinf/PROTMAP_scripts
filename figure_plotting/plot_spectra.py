#!/usr/bin/python3

from pyteomics import mzml
import subprocess
import os
import matplotlib.pyplot as plt
from spectrum_utils import plot
from spectrum_utils import spectrum
import json


SMALL_SIZE = 17
MEDIUM_SIZE = 25

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels


def load_data():
    selection_file = "nov_psm6"
    with open("./parameters.json", "r") as file_handle:
        parameters = json.load(file_handle)
        data_dir = parameters['data_dir']
        publication_dir = parameters["publication_dir"]

    accu_data_dir = data_dir + "/accumulated_data"
    cand_dir = data_dir + "/candidates"
    mzml_dir = data_dir + "/ms/SIHUMI/mzML/"

    with open(cand_dir + "/" + selection_file + "_list.json", "r") as file_handler:
        cand_list = json.load(file_handler)

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    return mzml_dir, cand_list, prot_dic, publication_dir


def plot_spectra(mzml_id, peptide, scan_id, mzml_dir, spec_pic_dir, psm_id):
    mzml_file = str(subprocess.check_output("find {} -name {}.mzML".format(mzml_dir, mzml_id), shell=True))
    mzml_file = mzml_file.replace("b'", "").replace("\\n'", "")
    with mzml.read(mzml_file) as reader:
        # auxiliary.print_tree(next(reader))
        for scan in reader:
            if not scan["index"] == int(scan_id) - 1:
                continue
            if "precursorList" not in scan.keys():
                print("no precursor list")
                return
            mz = scan['m/z array']
            intensity = scan['intensity array']
            identifier = scan['index']
            retention_time = float(scan['scanList']['scan'][0]["scan start time"]) * 60.0
            precursor_mz = scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
            precursor_charge = int(scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["charge state"])
            spec = spectrum.MsmsSpectrum(identifier, precursor_mz,
                                         precursor_charge, mz, intensity,
                                         retention_time=retention_time,
                                         peptide=peptide)
            min_mz, max_mz = 100, 1400
            fragment_tol_mass, fragment_tol_mode = 10, 'ppm'
            min_intensity, max_num_peaks = 0.05, 150
            scaling = 'root'
            ion_types = 'aby'
            spec = spec.set_mz_range(min_mz, max_mz)
            spec = spec.remove_precursor_peak(fragment_tol_mass,
                                              fragment_tol_mode)
            spec = spec.filter_intensity(min_intensity, max_num_peaks)
            spec = spec.scale_intensity(scaling)
            # spec = spec.annotate_peaks(fragment_tol_mass, fragment_tol_mode, ion_types)
            spec = spec.annotate_peptide_fragments(fragment_tol_mass,
                                                   fragment_tol_mode,
                                                   ion_types)
            plt.figure()
            plot.spectrum(spec, grid=False)
            mzml_id = os.path.splitext(os.path.split(mzml_file)[1])[0]
            plt.savefig("{}/{}_{}.svg".format(spec_pic_dir, mzml_id, psm_id), bbox_inches='tight')
            plt.close()
            print("print")
            return
        else:
            print("Scan not found")


def get_specs(candidate, spectra_id_list, prot_dic, spec_pic_dir, mzml_dir):
    if not os.path.exists(spec_pic_dir):
        os.makedirs(spec_pic_dir)

    for spectra_id in spectra_id_list:
        for psm in prot_dic['SIHUMI']['6frame'][candidate]:
            scan_id = psm["scan"]
            num = psm["num"]
            # uncomment if to get all spectra of candidate
            if scan_id != str(spectra_id[0]) or num != spectra_id[1]:
                continue
            psm_id = scan_id + "_" + str(num)
            mzml_id = psm["experiment"]
            peptide = psm["pep"]
            psm_id = scan_id + "_" + str(num)
            info_file = spec_pic_dir + mzml_id + "_" + psm_id + "_info.txt"
            with open(info_file, "w") as file_handle:
                file_handle.write("FDR: " + str(psm["FDR"]) + "%\n")
                file_handle.write("e-value: " + str(psm["e-value"]) + "\n")
                file_handle.write("X-Corr: " + str(psm["xcorr"]) + "\n")
                file_handle.write("id: " + psm_id + "\n")
                file_handle.write("experiment: " + psm["experiment"] + "\n")
                file_handle.write("Peptide: " + psm["pep"] + "\n")
            plot_spectra(mzml_id, peptide, scan_id, mzml_dir, spec_pic_dir, psm_id)


def main():
    mzml_dir, cand_list, prot_dic, publication_dir = load_data()
    spec_pic_root = publication_dir + "/figs/"

    candidate = "blautia|CP039126.1|146236|1|3946668-3946884"
    # selection of the specific spectra is not ambiguous but in this case clear
    spectra_id_list = [[3659, 1], [5175, 1], [7763, 1]]
    candidate_name = "nov_57"
    spec_pic_dir = spec_pic_root + candidate_name + "_genome_view/"
    get_specs(candidate, spectra_id_list, prot_dic,
              spec_pic_dir, mzml_dir)

    candidate = "blautia|CP039126.1|383626|-1|1942316-1942694"
    spectra_id_list = [[8045, 1], [8500, 1], [37240, 1]]
    candidate_name = "nov_5"
    spec_pic_dir = spec_pic_root + candidate_name + "_genome_view/"
    get_specs(candidate, spectra_id_list, prot_dic,
              spec_pic_dir, mzml_dir)


if __name__ == "__main__":
    main()


'''
mzml_dir, database_dic, cand_list, prot_dic = load_data()
spec_pic_root = "../figs/"

candidate = "blautia|CP039126.1|146236|1|3946668-3946884"
# selection of the specific spectra
spectra_id_list = [[3659, 1]]
candidate_name = "nov_57"
spec_pic_dir = spec_pic_root + candidate_name + "_genome_view/"
get_specs(candidate, spectra_id_list, prot_dic, spec_pic_dir, mzml_dir)

candidate = "blautia|CP039126.1|383626|-1|1942316-1942694"
spectra_id_list = [[8045, 1], [8500, 1], [42759, 1], [16840, 1], [18762, 1],
                   [37240, 1], [30384, 1], [12601, 1]]
candidate_name = "nov_5"
spec_pic_dir = "/homes/biertruck/john/public_html/"
'''
