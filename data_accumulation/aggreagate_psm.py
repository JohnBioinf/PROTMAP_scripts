#!/usr/bin/python3

from os import listdir, makedirs, path
from bisect import bisect_left
import json


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
        if all([decoy_fix in protein for protein in protein_list]):
            psm["decoy"] = True
        else:
            psm["decoy"] = False

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


def collect_results(ms_dir, crap_headers):
    print("loading results")
    # collect results and calculate FDR
    data_base_types = ["proteom", "6frame"]
    species_types = ["ecoli", "SIHUMI"]
    psm_base_dir_temp = ms_dir + "/{0}/psm/{1}/"

    database_dic = {}
    for species_type in species_types:
        print(species_type)
        database_dic[species_type] = {}
        for data_base in data_base_types:
            psm_base_dir = psm_base_dir_temp.format(species_type, data_base)
            print(data_base)
            database_dic[species_type][data_base] = {}
            for psm_exp in listdir(psm_base_dir):
                if not psm_exp.endswith(".txt"):
                    continue

                psm_dir = psm_base_dir + psm_exp
                experiment = psm_exp.replace(".txt", "")
                print(experiment)
                psms = load_psm(psm_dir, crap_headers, exp=experiment)
                database_dic[species_type][data_base][experiment] = psms
    return database_dic


def aggregate_prot(FDR_cut_off, database_dic):
    prot_dic = {}
    for species_type, species_database_dic in database_dic.items():
        print(species_type)
        prot_dic[species_type] = {}
        for data_base, DB_database_dic in species_database_dic.items():
            print(data_base)
            prot_dic[species_type][data_base] = {}
            for experiment, experiment_dic in DB_database_dic.items():
                print(experiment)
                for scan, psms in experiment_dic.items():
                    for psm in psms:
                        if psm["FDR"] <= FDR_cut_off:
                            for protein in psm["proteins"].split(","):
                                if "DECOY_" in protein:
                                    continue
                                if protein in prot_dic[species_type][data_base]:
                                    prot_dic[species_type][data_base][protein].append(psm)
                                else:
                                    prot_dic[species_type][data_base][protein] = [psm]

    return prot_dic


def main():
    print("loading data")
    FDR_cut_off = 1

    with open("./parameters.json", "r") as file_handle:
        parameters = json.load(file_handle)
        data_dir = parameters['data_dir']
        publication_dir = parameters['publication_dir']

    db_dir = data_dir + "/dbs"
    ms_dir = data_dir + "/ms"
    accu_data_dir = data_dir + "/accumulated_data"

    with open(db_dir + "/crap.fasta", "r") as file_handle:
        crap_headers = [l[1:-1] for l in file_handle if ">" in l]

    if not path.exists(accu_data_dir):
        makedirs(accu_data_dir)

    database_dic = collect_results(ms_dir, crap_headers)
    with open(accu_data_dir + "/database_dic.json", "w") as file_handle:
        json.dump(database_dic, file_handle)

    with open(accu_data_dir + "/database_dic.json", "r") as file_handle:
        database_dic = json.load(file_handle)

    prot_dic = aggregate_prot(FDR_cut_off, database_dic)
    with open(accu_data_dir + "/prot_dic.json", "w") as file_handle:
        json.dump(prot_dic, file_handle)

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    with open(publication_dir + "/ms_experiment_info.txt", "w") as file_handle:
        file_handle.write("Number of total experiments: " +
                          str(len(database_dic["ecoli"]["6frame"].keys()) +
                              len(database_dic["SIHUMI"]["6frame"].keys())) + "\n")
        stand_exp = len([1 for exp in database_dic["SIHUMI"]["6frame"].keys()
                         if "SIHUMI_standard" in exp])
        small_exp = len([1 for exp in database_dic["SIHUMI"]["6frame"].keys()
                         if "SIHUMI_small" in exp])
        file_handle.write("Number of SIHUMI expriments: " +
                          str(stand_exp + small_exp) + "\n")
        file_handle.write("Number of standard expriments: " + str(stand_exp) + "\n")
        file_handle.write("Number of small expriments: " + str(small_exp) + "\n")
        num_psm = 0
        for protein, psms in prot_dic["SIHUMI"]["6frame"].items():
            num_psm += len([psm for psm in psms if not psm["decoy"]])
        file_handle.write("Number of PSMs in SIHUMI with 6frame experiments: " +
                          str(num_psm) + "\n")

        num_psm = 0
        for protein, psms in prot_dic["ecoli"]["6frame"].items():
            num_psm += len([psm for psm in psms if not psm["decoy"]])
        for protein, psms in prot_dic["ecoli"]["proteom"].items():
            num_psm += len([psm for psm in psms if not psm["decoy"]])

        file_handle.write("Number of PSMs in ecoli data set: " + str(num_psm) + "\n")


if __name__ == "__main__":
    main()


'''
FDR_cut_off = 1

with open("./parameters.json", "r") as file_handle:
    data_dir = json.load(file_handle)['data_dir']

db_dir = data_dir + "/dbs"
ms_dir = data_dir + "/ms"
accu_data_dir = data_dir + "/accumulated_data"

with open(accu_data_dir + "/database_dic.json", "r") as file_handle:
    database_dic = json.load(file_handle)

with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
    prot_dic = json.load(file_handle)

print("Number of total experiments: " +
      str(len(database_dic["ecoli"]["6frame"].keys()) +
          len(database_dic["SIHUMI"]["6frame"].keys())))
stand_exp = len([1 for exp in database_dic["SIHUMI"]["6frame"].keys()
                 if "SIHUMI_standard" in exp])
small_exp = len([1 for exp in database_dic["SIHUMI"]["6frame"].keys()
                 if "SIHUMI_small" in exp])
print("Number of standard expriments: " + str(stand_exp))
print("Number of small expriments: " + str(small_exp))
num_psm = 0
for protein, psms in prot_dic["SIHUMI"]["6frame"].items():
    num_psm += len(psms)

print("Number of PSMs in SIHUMI experiments: " + str(num_psm))
for protein, psms in prot_dic["ecoli"]["6frame"].items():
    num_psm += len(psms)

print("Number of total PSMs: " + str(num_psm))
print("Number of analysed Spectra")
'''
