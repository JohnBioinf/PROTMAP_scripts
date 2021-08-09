#!/usr/bin/python3

'''Call like:
    ./pepG_v1.2.py /scratchsan2/john/MS-Spektren/psm/6frame_new_SIHUMI/SIHUMI/steffi/180228_StS_flowrate_Bd13_ctrl_rep1_id.txt /scratchsan2/john/MS-Spektren/res/6frame_new_sequenced_SIHUMI/SIHUMI/steffi/ 6frame
'''

import os
import re
import argparse
from argparse import RawTextHelpFormatter
from bisect import bisect_left
import json
from Bio import SeqIO


def scien_to_real(scien_num):
    if "E" in scien_num:
        return(float(scien_num.replace("E", "e")))
    else:
        return(float(scien_num))


def calc_FDR(decoy_psm, real_psm, e_cut_off):
    decoy_psm.sort()
    real_psm.sort()
    real_psm.append(e_cut_off)
    FDR_tuple = []
    for num_real_above, real_e in enumerate(real_psm):
        num_decoy_above = bisect_left(decoy_psm, real_e)
        FDR = round(float(num_decoy_above) / float(num_real_above + 1) * 100.0,
                    2)
        FDR_tuple.append((real_e, FDR))
    FDR_tuple = list(set(FDR_tuple))
    FDR_tuple.sort()
    return FDR_tuple


def load_psm(psm_dir, crap_headers, decoy_fix="DECOY_", e_cut_off=1,
             exp="unknown"):
    """input:
    path to file
    output:
    dic[scan_id] = [{psm}, {psm}, ...]
    keys of {psm}
    0   experiment
    1   scan
    2   e-value
    3   pep
    4   proteins
    5   xcorr
    6   delta_cn
    7   sp_score
    8   id
    9   num
    10  FDR
    """
    scan_pep_map = {}
    decoy_psm = []
    real_psm = []
    file_handle = open(psm_dir, "r")
    for line in file_handle:
        if line == "":
            continue
        if "CometVersion" in line:
            continue
        if "exp_neutral_mass" in line:
            continue
        line = line.split("\t")
        e_value = scien_to_real(line[5])
        protein_list = line[15].split(",")
        if any(crap in protein_list for crap in crap_headers):
            continue
        if e_value > e_cut_off:
            continue
        # check if all or none are decoy
        if len(set([decoy_fix in protein for protein in protein_list])) == 1:
            if decoy_fix in protein_list[0]:
                decoy_psm.append(e_value)
            else:
                real_psm.append(e_value)
        if all([decoy_fix in protein for protein in protein_list]):
            continue
        scan = line[0]
        psm = {}
        psm["experiment"] = exp
        psm["scan"] = scan
        psm["e-value"] = e_value
        psm["pep"] = line[11]
        psm["proteins"] = line[15]
        psm["xcorr"] = line[6]
        psm["delta_cn"] = line[7]
        psm["sp_score"] = line[8]
        psm["id"] = scan + "_" + line[1]
        psm["num"] = int(line[1])
        psm["FDR"] = None

        if scan in scan_pep_map:
            scan_pep_map[scan].append(psm)
        else:
            scan_pep_map[scan] = [psm]

    file_handle.close()
    FDR_tuple = calc_FDR(decoy_psm, real_psm, e_cut_off)
    e_value_list, FDR_list = zip(*FDR_tuple)
    for scan, psms in scan_pep_map.items():
        for psm in psms:
            e_value = psm["e-value"]
            FDR = FDR_list[bisect_left(e_value_list, e_value)]
            psm["FDR"] = FDR
    return scan_pep_map


def pep_hit_to_bed_line(protein, pep, DB, DB_name, experiment):
    blockCount = 1
    blockStarts = 0
    expIds = 1
    expCount = 1

    frame_color_dic = {}
    # frame 1 red
    frame_color_dic[1] = ["255,128,128", "255,26,26"]
    # frame 2 blue
    frame_color_dic[2] = ["128,128,255", "26,26,255"]
    # frame 3 green
    frame_color_dic[3] = ["153,230,153", "45,185,45"]
    # frame 4 purple
    frame_color_dic[4] = ["223,128,255", "172,0,230"]
    # frame 5 orange
    frame_color_dic[5] = ["255,191,128", "230,115,0"]
    # frame 6 turquoise
    frame_color_dic[6] = ["128,255,229", "0,230,184"]

    pep_seq = pep["pep"]
    FDR = pep["FDR"]
    psm_id = pep["id"]
    score = min(int(FDR * 10), 1000)
    num = pep["num"]

    protein_seq = DB[protein].seq
    chromosom = protein.split('|')[1]
    strand = protein.split("|")[3]
    org_start = int(protein.split("|")[4].split("-")[0])
    org_stop = int(protein.split("|")[4].split("-")[1])

    if (len(pep["proteins"].split(",")) == 1) and (num == 1):
        quality = 1
    else:
        quality = 0

    # MS cant discriminate between Isoleucin adn Leucin
    as_regex = re.compile("[IL]")
    pep_seq_regex = as_regex.sub("[IL]", pep_seq)
    pep_starts = [m.start() for m in re.finditer(pep_seq_regex,
                                                 str(protein_seq))]
    if(pep_starts == []):
        print("Peptide not in protein")
        print(pep)
        print(pep_seq)
        raise
    lines = []
    for pep_start in pep_starts:
        name = ("spec|" + DB_name + "|" + chromosom + "|" + experiment + "|" +
                psm_id)
        if strand == "-1":
            # start of protein in genome
            thickEnd = chromEnd = org_stop - pep_start * 3
            thickStart = chromStart = chromEnd - len(pep_seq) * 3
            frame = chromStart % 3 + 4
            strand = "-"
        else:
            # start of protein in genome
            thickStart = chromStart = org_start + pep_start * 3
            thickEnd = chromEnd = chromStart + len(pep_seq) * 3
            frame = chromStart % 3 + 1
            strand = "+"

        itemRgb = frame_color_dic[frame][quality]
        blockSizes = chromEnd - chromStart
        lines.append([chromosom, chromStart, chromEnd, name, score, strand,
                      thickStart, thickEnd, itemRgb, blockCount, blockSizes,
                      blockStarts, expCount, expIds])
    return lines


def psm_to_bed(scan_pep_map, FDR_cut_off, DB, experiment, args):
    track_template = ("track name=\"pepG {0}\" description="
                      "\"mapped peptides from {0} in {1}\" itemRgb=\"on\"")

    bed_dic = {}
    for scan, peptides in scan_pep_map.items():
        for pep in peptides:
            proteins = pep["proteins"].split(',')
            for protein in proteins:
                if "DECOY_" in protein:
                    continue
                if pep["FDR"] > FDR_cut_off:
                    continue
                bed_lines = pep_hit_to_bed_line(protein, pep, DB,
                                                args.protein_DB, experiment)
                species = protein.split("|")[0]
                if species in bed_dic.keys():
                    bed_dic[species] += bed_lines
                else:
                    track = track_template.format(experiment, species)
                    bed_dic[species] = [track]
                    bed_dic[species] += bed_lines
    return bed_dic


def print_results(bed_dic, experiment, args):
    for species in bed_dic.keys():
        output_dir = args.output_dir + "/" + species + "/"
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        output_file = output_dir + experiment + ".bed"
        if species in bed_dic:
            bed_lines = bed_dic[species]
            with open(output_file, "w") as file_handler:
                file_handler.write(bed_lines[0] + "\n")
                del bed_lines[0]
                bed_lines.sort(key=lambda x: (x[0], x[1]))
                file_handler.write("\n".join([" ".join(map(str, l)) for l in
                                              bed_lines]) + "\n")
        else:
            with open(output_file, "w") as file_handler:
                file_handler.write("")


def main():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)

    parser.add_argument("pep", help="comet result txt", type=str)
    parser.add_argument("output_dir", help="output directory", type=str)
    parser.add_argument("protein_DB", help=("the name of the database. Currently "
                                            "only 6frame or proteom."),
                        type=str, choices=["6frame", "proteom"])

    args = parser.parse_args()

    scan_pep_map = {}
    e_cut_off = 0.1
    FDR_cut_off = 1

    experiment = os.path.basename(args.pep).replace(".txt", "")
    print(experiment)
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    with open(data_dir + "/dbs/crap.fasta", "r") as file_handle:
        crap_headers = [l[1:-1] for l in file_handle if ">" in l]

    path_db = data_dir + "/dbs/SIHUMI_" + args.protein_DB + ".fasta"
    DB = SeqIO.index(path_db, "fasta")

    scan_pep_map = load_psm(args.pep, crap_headers, e_cut_off=e_cut_off)

    bed_dic = psm_to_bed(scan_pep_map, FDR_cut_off, DB, experiment, args)

    print_results(bed_dic, experiment, args)


if __name__ == "__main__":
    main()
