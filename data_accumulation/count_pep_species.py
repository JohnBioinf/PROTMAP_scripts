#!/usr/bin/python3

import json
from Bio import SeqIO


def load_data():
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    with open("./SIHUMI_info_dic.json", "r") as file_handle:
        SIHUMI_info_dic = json.load(file_handle)

    db_dir = data_dir + "/dbs"
    frame_dic = SeqIO.index(db_dir + "/SIHUMI_6frame.fasta", "fasta")

    accu_data_dir = data_dir + "/accumulated_data"

    with open(accu_data_dir + "/database_dic.json", "r") as file_handle:
        database_dic = json.load(file_handle)
    return database_dic, frame_dic, SIHUMI_info_dic, accu_data_dir


def count_spec_4_species(database_dic, FDR_cut_off):
    spec_count = {}
    min_pep_size = float("inf")
    for database_type, db_dic in database_dic["SIHUMI"].items():
        print(database_type)
        for experiment, scan_dic in db_dic.items():
            for scan, psms in scan_dic.items():
                for psm in psms:
                    if psm["decoy"]:
                        continue
                    if psm["FDR"] > FDR_cut_off:
                        continue
                    if "DECOY_" in psm["proteins"]:
                        continue
                    if len(psm["proteins"].split(",")) > 1:
                        continue
                    if len(psm["pep"]) < min_pep_size:
                        min_pep_size = len(psm["pep"])
                    species = psm["proteins"].split("|")[0]
                    if species in spec_count:
                        spec_count[species] += 1
                    else:
                        spec_count[species] = 1

    return spec_count, min_pep_size


def get_size_of_proteogenome(frame_dic, SIHUMI_info_dic, min_pep_size):
    size_proteogenome = {}
    for spec in SIHUMI_info_dic.keys():
        sum = 0
        for orf, prot_seq in frame_dic.items():
            # if orf is smaller than possible pep not possible to map
            # further likelihood of size increases with size
            size = len(prot_seq.seq) - min_pep_size
            if size > 0:
                sum += size
        size_proteogenome[spec] = sum

    return size_proteogenome


def main():
    FDR_cut_off = 1
    database_dic, frame_dic, SIHUMI_info_dic, accu_data_dir = load_data()

    spec_count, min_pep_size = count_spec_4_species(database_dic, FDR_cut_off)

    size_proteogenome = get_size_of_proteogenome(frame_dic, SIHUMI_info_dic,
                                                 min_pep_size)

    # normalize
    res = ""
    for species, count in spec_count.items():
        count_norm = count / size_proteogenome[species]
        res += "{}\t{}\n".format(SIHUMI_info_dic[species]["short_name"],
                                 count_norm)

    with open(accu_data_dir + "/pep_species_ratio.tsv", "w") as file_handler:
        file_handler.write(res)


if __name__ == "__main__":
    main()


'''
FDR_cut_off = 1
database_dic, frame_dic, SIHUMI_info_dic, accu_data_dir = load_data()

spec_count, min_pep_size = count_spec_4_species(database_dic, FDR_cut_off)

size_proteogenome = get_size_of_protegenom(frame_dic, SIHUMI_info_dic,
                                           min_pep_size)

res = ""
for species, count in spec_count.items():
    count_norm = count / size_proteogenome[species]
    res += "{}\t{}\n".format(SIHUMI_info_dic[species]["short_name"],
                             count_norm)

with open(accu_data_dir + "/pep_species_ratio.tsv", "w") as file_handler:
    file_handler.write(res)
'''
