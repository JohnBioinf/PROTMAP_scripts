#!/usr/bin/python3

import json
from itertools import chain


def load_data():
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    db_dir = data_dir + "/dbs"
    accu_data_dir = data_dir + "/accumulated_data"

    with open(db_dir + "/CDS_to_orf.json", "r") as file_handle:
        CDS_to_orf = json.load(file_handle)

    with open(db_dir + "/orf_to_CDS.json", "r") as file_handle:
        orf_to_CDS = json.load(file_handle)

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    return accu_data_dir, orf_to_CDS, CDS_to_orf, prot_dic


def data_4_venn(orf_to_CDS, CDS_to_orf, accu_data_dir, prot_dic):
    with open(accu_data_dir + "/all_CDS", "w") as file_handler:
        all_CDS = [orf[0] for cds, orf in CDS_to_orf["ecoli"].items()]
        all_CDS = list(set(list(chain.from_iterable(all_CDS))))
        file_handler.write("\n".join(all_CDS))

    ks = [1, 6, 10]
    species_type = "ecoli"
    for k in ks:
        with open(accu_data_dir + "/k_{}_6frame".format(k), "w") as file_handler:
            the_set = [prot
                       for prot, psm in prot_dic[species_type]["6frame"].items()
                       if len(psm) >= k]
            file_handler.write("\n".join(the_set))
        with open(accu_data_dir + "/k_{}_proteom".format(k), "w") as file_handler:
            the_set = [CDS_to_orf[prot.split("|")[0]][prot][0]
                       for prot, psm in prot_dic[species_type]["proteom"].items()
                       if len(psm) >= k]
            the_set = list(chain.from_iterable(the_set))
            file_handler.write("\n".join(the_set))

    for k in ks:
        orf_proteom = [CDS_to_orf[prot.split("|")[0]][prot][0]
                       for prot, psm in prot_dic[species_type]["proteom"].items()
                       if len(psm) >= k]
        orf_proteom = set(chain.from_iterable(orf_proteom))
        orf_6frame = set([prot
                          for prot, psm in prot_dic[species_type]["6frame"].items()
                          if len(psm) >= k])
        only_6frame = orf_6frame - orf_proteom
        only_proteom = [orf_to_CDS[species_type][orf][0]
                        for orf in orf_proteom - orf_6frame]
        only_6frame_noCDS = [prot for prot in only_6frame
                             if orf_to_CDS["ecoli"][prot][0] is None]

        with open(accu_data_dir + "/k_{}_scores_6frame.txt".format(k), "w") as file_handler:
            scores_6frame = []
            for prot, psms in prot_dic[species_type]["6frame"].items():
                if len(psms) < k:
                    continue
                scores_6frame.append("\n".join([str(psm["e-value"]) for psm in psms]))
            file_handler.write("\n".join(scores_6frame))

        with open(accu_data_dir + "/k_{}_scores_proteom.txt".format(k), "w") as file_handler:
            scores_proteom = []
            for prot, psms in prot_dic[species_type]["proteom"].items():
                if len(psms) < k:
                    continue
                scores_proteom.append("\n".join([str(psm["e-value"]) for psm in psms]))
            file_handler.write("\n".join(scores_proteom))

        with open(accu_data_dir + "/k_{}_scores_only_proteom.txt".format(k), "w") as file_handler:
            scores_only_proteom = []
            for prot, psms in prot_dic[species_type]["proteom"].items():
                if (len(psms) >= k) and (prot in only_proteom):
                    scores_only_proteom.append("\n".join([str(psm["e-value"]) for psm in
                                                          psms]))
            file_handler.write("\n".join(scores_only_proteom))

        with open(accu_data_dir + "/k_{}_scores_only_6frame.txt".format(k), "w") as file_handler:
            scores_only_6frame = []
            for prot, psms in prot_dic[species_type]["6frame"].items():
                if (len(psms) >= k) and (prot in only_6frame):
                    scores_only_6frame.append("\n".join([str(psm["e-value"]) for psm in
                                                         psms]))
            file_handler.write("\n".join(scores_only_6frame))

        with open(accu_data_dir + "/k_{}_scores_only_6frame_noCDS.txt".format(k), "w") as file_handler:
            scores = []
            for prot, psms in prot_dic[species_type]["6frame"].items():
                if (len(psms) >= k) and (prot in only_6frame_noCDS):
                    scores.append("\n".join([str(psm["e-value"]) for psm in psms]))
            file_handler.write("\n".join(scores))


def data_4_edgecases(prot_dic, accu_data_dir):
    species_type = "ecoli"
    with open("./not_annotated_k10_ecoli.csv", "r") as file_handle:
        res = []
        for line in file_handle:
            candidate = line.split(",")[2]
            category = line.split(",")[4]
            e_vals = [str(psm["e-value"]) for psm in prot_dic[species_type]["6frame"][candidate]]
            res += [[e, category] for e in e_vals]

    with open(accu_data_dir + "/scores_not_annotated_k10_ecoli.csv", "w") as file_handle:
        file_handle.write("\n".join([",".join(line) for line in res]))

    # with open("./early_starts.csv", "r") as file_handle:
    #     res = []
    #     for line in file_handle:
    #         candidate = line.split(",")[2]
    #         category = line.split(",")[4]
    #         e_vals = [str(psm["e-value"]) for psm in prot_dic[species_type]["6frame"][candidate]]
    #         res += [[e, category] for e in e_vals]


def main():
    accu_data_dir, orf_to_CDS, CDS_to_orf, prot_dic = load_data()

    data_4_venn(orf_to_CDS, CDS_to_orf, accu_data_dir, prot_dic)

    data_4_edgecases(prot_dic, accu_data_dir)


if __name__ == "__main__":
    main()


'''
accu_data_dir, orf_to_CDS, CDS_to_orf, prot_dic = load_data()
k = 1
species_type = "ecoli"
orf_proteom = [CDS_to_orf[prot.split("|")[0]][prot][0]
               for prot, psm in prot_dic[species_type]["proteom"].items()
               if len(psm) >= k]
orf_proteom = set(chain.from_iterable(orf_proteom))
orf_6frame = set([prot
                  for prot, psm in prot_dic[species_type]["6frame"].items()
                  if len(psm) >= k])
only_6frame = orf_6frame - orf_proteom
only_proteom = [orf_to_CDS[species_type][orf][0]
                for orf in orf_proteom - orf_6frame]
scores_only_proteom = []
for prot, psms in prot_dic[species_type]["proteom"].items():
    if (len(psms) >= k) and (prot in only_proteom):
        scores_only_proteom.append("\n".join([str(psm["e-value"]) for psm in
                                              psms]))
'''
