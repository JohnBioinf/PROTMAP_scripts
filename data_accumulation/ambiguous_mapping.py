#!/usr/bin/python3

import json


def load_data():
    with open("./parameters.json", "r") as file_handle:
        parameters = json.load(file_handle)
        data_dir = parameters['data_dir']
        publication_dir = parameters["publication_dir"]

    accu_data_dir = data_dir + "/accumulated_data"

    with open(accu_data_dir + "/database_dic.json", "r") as file_handle:
        database_dic = json.load(file_handle)

    return database_dic, publication_dir


def make_ambiguous_dic(database_dic, FDR_cut_off):
    ambiguous_dic = {}
    for database_type, db_dic in database_dic["SIHUMI"].items():
        print(database_type)
        for experiment, scan_dic in db_dic.items():
            for scan, psms in scan_dic.items():
                for psm in psms:
                    if "DECOY_" in psm["proteins"]:
                        continue
                    key = database_type + "_" + str(len(psm["proteins"].split(",")))
                    if key in ambiguous_dic:
                        ambiguous_dic[key] += 1
                    else:
                        ambiguous_dic[key] = 1

    return ambiguous_dic


def make_table(ambiguous_dic, max_num, sum_proteom, sum_6frame):
    table_tex_list = [None] * max_num
    for i in range(1, max_num + 1):
        num_prot = str(round((ambiguous_dic.get("proteom_" + str(i), 0) /
                              sum_proteom) * 100, 4)) + " \\%"
        num_6frame = str(round((ambiguous_dic.get("6frame_" + str(i), 0) /
                                sum_6frame) * 100, 4)) + " \\%"
        table_tex_list[i - 1] = map(str, [i, num_prot, num_6frame])

    return " \\\\\n".join([" & ".join(e) for e in table_tex_list]) + "\n"


def main():
    database_dic, publication_dir = load_data()

    ambiguous_dic = make_ambiguous_dic(database_dic)

    max_num = max([int(k.split("_")[1]) for k in ambiguous_dic.keys()])
    sum_proteom = sum([int(v) for k, v in ambiguous_dic.items()
                       if "proteom" in k])
    sum_6frame = sum([int(v) for k, v in ambiguous_dic.items()
                      if "6frame" in k])

    table_tex = make_table(ambiguous_dic, max_num, sum_proteom, sum_6frame)

    with open(publication_dir + "/ambiguous_mapping.tex", "w") as file_handler:
        file_handler.write(table_tex)


if __name__ == "__main__":
    main()


'''
FDR_cut_off = 1

database_dic, accu_data_dir = load_data()

ambiguous_dic = make_ambiguous_dic(database_dic, FDR_cut_off)

max_num = max([int(k.split("_")[1]) for k in ambiguous_dic.keys()])
sum_proteom = sum([int(v) for k, v in ambiguous_dic.items()
                   if "proteom" in k])
sum_6frame = sum([int(v) for k, v in ambiguous_dic.items()
                  if "6frame" in k])

table_tex = make_table(ambiguous_dic, max_num, sum_proteom, sum_6frame)
'''
