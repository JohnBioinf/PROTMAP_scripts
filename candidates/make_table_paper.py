#!/usr/bin/python3

import json
import sys
from Bio import SeqIO


def load_data(selection_file):
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    db_dir = data_dir + "/dbs"
    genome_dir = data_dir + "/genome"
    accu_data_dir = data_dir + "/accumulated_data"

    with open(db_dir + "/orf_to_CDS.json", "r") as file_handle:
        orf_to_CDS = json.load(file_handle)

    with open(genome_dir + "/validation_dic.json", "r") as file_handle:
        validation_dic = json.load(file_handle)

    proteom_dic = SeqIO.index(db_dir + "/SIHUMI_6frame.fasta", "fasta")

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    return (proteom_dic, validation_dic, prot_dic, orf_to_CDS)


def collect_infos(prot_dic, species_list, orf_to_CDS, validation_dic):
    result_dic = {}
    k_up = 10
    k_low = 6

    for s in species_list:
        result_dic[s] = {}
        for k in [k_up, k_low]:
            for prot_type in ["nov", "hypo", "known"]:
                result_dic[s][prot_type + "_" + str(k)] = set()

    for protein, psms in prot_dic["SIHUMI"]["6frame"].items():
        species = protein.split("|")[0]
        if len(protein.split(";")) > 1:
            continue
        if orf_to_CDS[species][protein][0] is None:
            prot_type = "nov"
        else:
            protein_id = orf_to_CDS[species][protein][0].split("|")[2]
            validation_level = int(validation_dic[species][protein_id][0])
            if validation_level < 3:
                prot_type = "hypo"
            else:
                prot_type = "known"
        num_good_psms = len([1 for psm in psms if psm["FDR"] <= 1])
        if num_good_psms < k_low:
            continue
        k = k_up if num_good_psms >= k_up else k_low

        if k == k_up:
            key = prot_type + "_" + str(k_up)
            result_dic[species][key].add(protein)
        key = prot_type + "_" + str(k_low)
        result_dic[species][key].add(protein)

    return result_dic


def print_table(result_dic, species_list, proteom_dic):
    k_up = 10
    k_low = 6
    species_num_anno_dic = {}
    for protein in proteom_dic:
        species = protein.split("|")[0]
        if species not in species_num_anno_dic:
            species_num_anno_dic[species] = 0
        species_num_anno_dic[species] += 1

    species_science_name_dic = {"bact": "\\emph{B.\\ theta.}",
                                "blautia": "\\emph{B.\\ producta}",
                                "ecoli": "\\emph{E.\\ coli}",
                                "anaero": "\\emph{A.\\ caccae}",
                                "bifi": "\\emph{B.\\ longum}",
                                "ery": "\\emph{E.\\ ramosum}",
                                "lacto": "\\emph{L.\\ plantarum}",
                                "clostri": "\\emph{C.\\ butyricum}"}

    res = ""
    res += str(k_up) + "\n"
    for species in species_list:
        results = result_dic[species]
        tab_temp = "{} & {} & {} & {} & {}\\\\"
        res += tab_temp.format(species_science_name_dic[species],
                               len(results["nov_" + str(k_up)]),
                               len(results["hypo_" + str(k_up)]),
                               len(results["known_" + str(k_up)]),
                               round((len(results["known_" + str(k_up)]) +
                                      len(results["hypo_" + str(k_up)])) /
                                     species_num_anno_dic[species] * 100, 1))

    res += str(k_low) + "\n"
    for species in species_list:
        results = result_dic[species]
        tab_temp = "{} & {} & {} & {} & {}\\\\"
        res += tab_temp.format(species_science_name_dic[species],
                               len(results["nov_" + str(k_low)]),
                               len(results["hypo_" + str(k_low)]),
                               len(results["known_" + str(k_low)]),
                               round((len(results["known_" + str(k_low)]) +
                                      len(results["hypo_" + str(k_low)])) /
                                     species_num_anno_dic[species] * 100, 1))

    return res


def main():
    selection_cut_off = int(sys.argv[1])
    selection_file = "nov_psm" + str(selection_cut_off)
    # change order for order of table
    species_list = ['bact', 'blautia', 'ecoli', 'ery', 'bifi', 'anaero',
                    'lacto', 'clostri']
    (proteom_dic, validation_dic, prot_dic,
     orf_to_CDS) = load_data(selection_file)

    result_dic = collect_infos(prot_dic, species_list,  orf_to_CDS,
                               validation_dic)

    with open("../protein_detected.tex", "w") as file_handle:
        file_handle.write(print_table(result_dic, species_list, proteom_dic))


if __name__ == "__main__":
    main()

'''
selection_cut_off = 10
selection_file = "nov_psm" + str(selection_cut_off)

(proteom_dic, SIHUMI_info_dic, validation_dic, prot_dic,
 orf_to_CDS) = load_data(selection_file)

result_dic = collect_infos(prot_dic, SIHUMI_info_dic,  orf_to_CDS,
                           validation_dic)

print_table(result_dic, SIHUMI_info_dic)
'''
