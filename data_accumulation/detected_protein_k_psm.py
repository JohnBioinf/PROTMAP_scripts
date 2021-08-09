#!/usr/bin/python3

import json


def make_results(prot_dic, orf_to_CDS):
    res = ""
    for k in range(1, 21):
        orf_proteom = set([prot for prot, psm in prot_dic["ecoli"]["proteom"].items()
                           if len(psm) >= k])
        orf_6frame = set([prot for prot, psm in prot_dic["ecoli"]["6frame"].items()
                          if len(psm) >= k])
        CDS_orf_6frame = set()
        for orf in orf_6frame:
            species = orf.split("|")[0]
            CDS = (orf_to_CDS[species][orf][0]
                   if orf_to_CDS[species][orf][0] is not None else orf)
            CDS_orf_6frame.add(CDS)

        only_6frame = CDS_orf_6frame - orf_proteom
        only_proteom = orf_proteom - CDS_orf_6frame
        # intersect = orf_proteom.intersection(CDS_orf_6frame)
        res += "{}\t{}\t{}\n".format(k, len(only_6frame), "only 6frame ecoli")
        res += "{}\t{}\t{}\n".format(k, len(only_proteom), "only proteom ecoli")
        # res += "{}\t{}\t{}\n".format(k, len(intersect), "intersect ecoli")

    for k in range(1, 21):
        orf_proteom = set([prot for prot, psm in prot_dic["SIHUMI"]["proteom"].items()
                           if len(psm) >= k])
        orf_6frame = set([prot for prot, psm in prot_dic["SIHUMI"]["6frame"].items()
                          if len(psm) >= k])
        CDS_orf_6frame = set()
        for orf in orf_6frame:
            species = orf.split("|")[0]
            CDS = (orf_to_CDS[species][orf][0]
                   if orf_to_CDS[species][orf][0] is not None else orf)
            CDS_orf_6frame.add(CDS)

        only_6frame = CDS_orf_6frame - orf_proteom
        only_proteom = orf_proteom - CDS_orf_6frame
        # intersect = orf_proteom.intersection(CDS_orf_6frame)
        res += "{}\t{}\t{}\n".format(k, len(only_6frame), "only 6frame SIHUMI")
        res += "{}\t{}\t{}\n".format(k, len(only_proteom), "only proteom SIHUMI")
        # res += "{}\t{}\t{}\n".format(k, len(intersect), "intersect SIHUMI")

    return res


def load_data():
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    accu_data_dir = data_dir + "/accumulated_data"
    db_dir = data_dir + "/dbs"

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    with open(db_dir + "/orf_to_CDS.json", "r") as file_handle:
        orf_to_CDS = json.load(file_handle)

    return prot_dic, orf_to_CDS, accu_data_dir


def main():
    prot_dic, orf_to_CDS, accu_data_dir = load_data()

    with open(accu_data_dir + "/detected_protein_k_psm.tsv", "w") as file_handle:
        file_handle.write(make_results(prot_dic, orf_to_CDS))


if __name__ == "__main__":
    main()

'''
prot_dic, orf_to_CDS = load_data()
make_results(prot_dic, orf_to_CDS)
'''
