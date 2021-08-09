#!/usr/bin/python3

import json
import csv
import sys


def load_data(selection_file):
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    cand_dir = data_dir + "/candidates"
    accu_data_dir = data_dir + "/accumulated_data"

    with open(cand_dir + "/" + selection_file + "_list.json", "r") as file_handle:
        candidate_list = json.load(file_handle)

    with open(cand_dir + "/" + selection_file + "_info_dic.json", "r") as file_handle:
        info_dic = json.load(file_handle)

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    return info_dic, candidate_list, prot_dic, accu_data_dir


def build_res(prot_dic, candidate_list, info_dic):
    res = []
    for species_type, species_prot_dic in prot_dic.items():
        print(species_type)
        for protein, psms in species_prot_dic["6frame"].items():
            if protein not in candidate_list:
                continue
            info = info_dic[protein]
            for psm in psms:
                if info["blast_category"] != "Novel" and info["blast_category"] != "Hypothetical":
                    blast_category = "Known"
                else:
                    blast_category = info["blast_category"]
                if info["start_codon"] == "ATG":
                    start_codon = "canonical"
                else:
                    start_codon = "not_canonical"

                res.append([psm["e-value"],
                            info["mean_top_psms"],
                            start_codon,
                            info["unique_psms"],
                            info["genomic_context_translated"],
                            info["genomic_context_category"],
                            info["genomic_context_strand"],
                            info["score_trans"],
                            blast_category,
                            info["name"],
                            info["candidate"]])

    return res


def main():
    selection_cut_off = int(sys.argv[1])
    selection_file = "nov_psm" + str(selection_cut_off)
    info_dic, candidate_list, prot_dic, accu_data_dir = load_data(selection_file)

    res = build_res(prot_dic, candidate_list, info_dic)

    with open(accu_data_dir + "/" + selection_file + "_scatterplot_s-hat.csv", "w") as file_handle:
        csv.writer(file_handle, delimiter=",").writerows(res)


if __name__ == "__main__":
    main()


'''
selection_file = "nov_psm6"
info_dic, candidate_list, prot_dic, accu_data_dir = load_data(selection_file)

res = build_res(prot_dic, candidate_list, info_dic)

with open(accu_data_dir + "/" + selection_file + "_scatterplot_mean.csv", "w") as file_handle:
    csv.writer(file_handle, delimiter=",").writerows(res)

'''
