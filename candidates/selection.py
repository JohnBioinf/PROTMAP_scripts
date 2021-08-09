#!/usr/bin/python3

import json
from Bio import SeqIO
import os
import sys


def selection(prot_dic, orf_to_CDS, selection_cut_off):
    cand_list = []
    for protein, psms in prot_dic["SIHUMI"]["6frame"].items():
        species = protein.split("|")[0]
        if orf_to_CDS[species][protein][1] is not None:
            continue
        num_good_psm = len([1 for psm in psms if psm["FDR"] <= 1])
        cand_list.append((num_good_psm, protein))

    cand_list = [e for e in cand_list if e[0] >= selection_cut_off]
    cand_list.sort(reverse=True)
    return [e[1] for e in cand_list]


def main():
    selection_cut_off = int(sys.argv[1])
    selection_file = "nov_psm" + str(selection_cut_off)

    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    db_dir = data_dir + "/dbs"
    with open(db_dir + "/orf_to_CDS.json", "r") as file_handle:
        orf_to_CDS = json.load(file_handle)

    genome_dir = data_dir + "/genome"
    genome = {}
    for genbank_file in os.listdir(genome_dir):
        if genbank_file.endswith(".gbk"):
            species = genbank_file.split(".")[0]
            genome_path = genome_dir + "/" + genbank_file
            genome[species] = SeqIO.to_dict(SeqIO.parse(genome_path, "genbank"))

    accu_data_dir = data_dir + "/accumulated_data"

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    cand_list = selection(prot_dic, orf_to_CDS, selection_cut_off)

    cand_data_dir = data_dir + "/candidates"
    with open(cand_data_dir + "/" + selection_file + "_list.json", "w") as file_handler:
        json.dump(cand_list, file_handler)


if __name__ == "__main__":
    main()


'''
'''
