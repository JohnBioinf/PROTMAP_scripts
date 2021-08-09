#!/usr/bin/python3

import json
from Bio import SeqIO


def load_data():
    with open("./parameters.json", "r") as file_handle:
        parameters = json.load(file_handle)
        data_dir = parameters['data_dir']
        publication_dir = parameters["publication_dir"]

    db_dir = data_dir + "/dbs"
    genome_dir = data_dir + "/genome"
    accu_data_dir = data_dir + "/accumulated_data"

    with open(db_dir + "/orf_to_CDS.json", "r") as file_handle:
        orf_to_CDS = json.load(file_handle)

    with open(genome_dir + "/validation_dic.json", "r") as file_handle:
        validation_dic = json.load(file_handle)

    proteom_dic = SeqIO.index(db_dir + "/SIHUMI_proteom.fasta", "fasta")

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    return proteom_dic, validation_dic, prot_dic, orf_to_CDS, publication_dir


def collect_infos(prot_dic, species_list, orf_to_CDS, validation_dic):
    result_dic = {}
    k_up = 10
    k_low = 6

    for s in species_list:
        result_dic[s] = {}
        for db_type in prot_dic["SIHUMI"]:
            for k in [k_up, k_low]:
                for prot_type in ["nov", "hypo", "known"]:
                    result_dic[s][prot_type + "_" + db_type + "_" +
                                  str(k)] = set()

    for db_type, db_prot_dic in prot_dic["SIHUMI"].items():
        for protein, psms in db_prot_dic.items():
            species = protein.split("|")[0]
            if len(protein.split(";")) > 1:
                continue
            prot_type = ""
            if db_type == "proteom":
                protein_id = protein.split("|")[2]
            else:
                if orf_to_CDS[species][protein][0] is None:
                    prot_type = "nov"
                else:
                    protein_id = orf_to_CDS[species][protein][0].split("|")[2]

            if prot_type != "nov":
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
                key = prot_type + "_" + db_type + "_" + str(k_up)
                result_dic[species][key].add(protein)
            key = prot_type + "_" + db_type + "_" + str(k_low)
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

    # main table
    tab_temp = "{} & {} & {} & {} & {}\\\\\n"
    res = ""
    for k in [k_up, k_low]:
        res += str(k) + "\n"
        for species in species_list:
            results = result_dic[species]
            res += tab_temp.format(species_science_name_dic[species],
                                   len(results["nov_6frame_" + str(k)]),
                                   len(results["hypo_6frame_" + str(k)]),
                                   len(results["known_6frame_" + str(k)]),
                                   round((len(results["known_6frame_" + str(k)]) +
                                          len(results["hypo_6frame_" + str(k)])) /
                                         species_num_anno_dic[species] * 100, 1))

    # supplement table
    tab_temp = "{} & {} & {} & {} & {} & {} & {} & {}\\\\\n"
    for k in [k_up, k_low]:
        res += str(k) + "\n"
        for species in species_list:
            results = result_dic[species]
            sum_6frame = sum([len(v) for k, v in results.items()
                              if "6frame" in k and "_" + str(k_low) in k])
            sum_proteom = sum([len(v) for k, v in results.items()
                              if "proteom" in k and "_" + str(k_low) in k])
            res += tab_temp.format(species_science_name_dic[species],
                                   len(results["nov_6frame_" + str(k)]),
                                   len(results["hypo_6frame_" + str(k)]),
                                   len(results["hypo_proteom_" + str(k)]),
                                   len(results["known_6frame_" + str(k)]),
                                   len(results["known_proteom_" + str(k)]),
                                   sum_proteom,
                                   sum_6frame)
    return res


def main():
    # change order for order of table
    species_list = ['bact', 'blautia', 'ecoli', 'ery', 'bifi', 'anaero',
                    'lacto', 'clostri']
    (proteom_dic, validation_dic, prot_dic,
     orf_to_CDS, publication_dir) = load_data()

    result_dic = collect_infos(prot_dic, species_list,  orf_to_CDS,
                               validation_dic)

    with open(publication_dir + "/protein_detected.tex", "w") as file_handle:
        file_handle.write(print_table(result_dic, species_list, proteom_dic))


if __name__ == "__main__":
    main()

'''
species_list = ['bact', 'blautia', 'ecoli', 'ery', 'bifi', 'anaero', 'lacto', 'clostri']
(proteom_dic, validation_dic, prot_dic, orf_to_CDS, CDS_to_orf) = load_data()

result_dic = collect_infos(prot_dic, species_list,  orf_to_CDS, validation_dic, CDS_to_orf)

print_table(result_dic, species_list, proteom_dic)
with open("../protein_detected.tex", "w") as file_handle:
    file_handle.write(print_table(result_dic, species_list, proteom_dic))
'''
