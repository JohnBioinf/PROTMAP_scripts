#!/usr/bin/python3
'''The two csv that are loaded where manual build'''

import json


def main():
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    accu_data_dir = data_dir + "/accumulated_data"

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    with open("./not_annotated_k10_ecoli.csv", "r") as file_handle:
        res = []
        for line in file_handle:
            candidate = line.split(",")[2]
            category = line.split(",")[4]
            e_vals = []
            if candidate in prot_dic["ecoli"]["6frame"]:
                e_vals += [str(psm["e-value"])
                           for psm in prot_dic["ecoli"]["6frame"][candidate]]
            res += [[e, category] for e in e_vals]

    with open(accu_data_dir + "/scores_not_annotated_k10_ecoli.csv", "w") as file_handle:
        file_handle.write("".join([",".join(line) for line in res]))

    with open("./early_starts.csv", "r") as file_handle:
        res = []
        for line in file_handle:
            candidate = line.split(",")[2]
            category = line.split(",")[4]
            if candidate in prot_dic["ecoli"]["6frame"]:
                e_vals += [str(psm["e-value"])
                           for psm in prot_dic["ecoli"]["6frame"][candidate]]
            if candidate in prot_dic["SIHUMI"]["6frame"]:
                e_vals += [str(psm["e-value"])
                           for psm in prot_dic["SIHUMI"]["6frame"][candidate]]
            res += [[e, category] for e in e_vals]

    with open(accu_data_dir + "/scores_early_starts.csv", "w") as file_handle:
        file_handle.write("".join([",".join(line) for line in res]))


if __name__ == "__main__":
    main()

'''
with open("./parameters.json", "r") as file_handle:
    data_dir = json.load(file_handle)['data_dir']

accu_data_dir = data_dir + "/accumulated_data"

with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
    prot_dic = json.load(file_handle)

with open("../not_annotated_k10_ecoli.csv", "r") as file_handle:
    res = []
    for line in file_handle:
        candidate = line.split(",")[2]
        category = line.split(",")[4]
        e_vals = []
        if candidate in prot_dic["ecoli"]["6frame"]:
            e_vals += [str(psm["e-value"])
                       for psm in prot_dic["ecoli"]["6frame"][candidate]]
        if candidate in prot_dic["SIHUMI"]["6frame"]:
            e_vals += [str(psm["e-value"])
                       for psm in prot_dic["SIHUMI"]["6frame"][candidate]]
        res += [[e, category] for e in e_vals]
'''
