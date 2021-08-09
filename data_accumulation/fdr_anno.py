#!/usr/bin/python3


from os import listdir, path
import json
import SIHUMI
from intervaltree import IntervalTree
from bisect import bisect_left
import re


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


def scien_to_real(scien_num):
    if "E" in scien_num:
        return(float(scien_num.replace("E", "e")))
    else:
        return(float(scien_num))


def load_psm(file_content, crap_headers, decoy_fix="DECOY_", e_cut_off=1,
             exp="unknown"):
    """input:
    str of file_content
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
    for line in file_content.split("\n"):
        if line == "":
            continue
        if "CometVersion" in line:
            continue
        if "exp_neutral_mass" in line:
            continue
        line = line.split("\t")
        e_value = scien_to_real(line[6])
        protein_list = line[16].split(",")
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
        if all([decoy_fix in protein for protein in protein_list]):
            continue
        scan = line[1]
        psm = {}
        psm["experiment"] = exp
        psm["scan"] = scan
        psm["e-value"] = e_value
        psm["pep"] = line[12]
        psm["proteins"] = line[16]
        psm["xcorr"] = line[7]
        psm["delta_cn"] = line[8]
        psm["sp_score"] = line[9]
        psm["id"] = line[0]
        psm["num"] = line[2]
        psm["FDR"] = None

        if scan in scan_pep_map:
            scan_pep_map[scan].append(psm)
        else:
            scan_pep_map[scan] = [psm]
    FDR_tuple = calc_FDR(decoy_psm, real_psm, e_cut_off)
    e_value_list, FDR_list = zip(*FDR_tuple)
    for scan, psms in scan_pep_map.items():
        for psm in psms:
            e_value = psm["e-value"]
            FDR = FDR_list[bisect_left(e_value_list, e_value)]
            psm["FDR"] = FDR
    return scan_pep_map


def collect_exp(psm_base_dir, crap_headers, e_cut_off):
    exp_dic = {}
    for psm_exp in listdir(psm_base_dir):
        psm_dir = psm_base_dir + psm_exp
        if not path.isdir(psm_dir):
            continue
        print(psm_dir)
        for psm in listdir(psm_dir):
            if psm.endswith("_id.txt"):
                experiment = psm.replace("_id.txt", "")
                if "Anne_Mendler" not in experiment:
                    continue
                print(experiment)
                file_content = ""
                with open(psm_dir + "/" + psm, 'r') as file_handler:
                    file_content += file_handler.read()
                exp_dic[experiment] = load_psm(file_content, crap_headers,
                                               e_cut_off=e_cut_off)
    return exp_dic


def bed2IntervalTree_dic(bed, IntervalTree_dic=None):
    """ input bed_file as list of strings for each line
    output dic with IntvalTree for each Chromosom/ Contig"""

    if IntervalTree_dic is None:
        IntervalTree_dic = {}

    for anno in bed:
        if ("track name" in anno) or (anno == ""):
            continue
        anno = anno.split(" ")
        if anno[3].split("|")[2] != "CDS":
            continue
        contig = anno[0]
        start = int(anno[1])
        stop = int(anno[2])
        name = anno[3]
        if contig not in IntervalTree_dic.keys():
            IntervalTree_dic[contig] = IntervalTree()
        # + 1 must be done as intevall are stored like so [start, stop)
        IntervalTree_dic[contig][start:stop + 1] = name
    return IntervalTree_dic


def collect_anno(scan_pep_map, IntervalTree_dic, orf_cds_dic, frame6,
                 value="e-value"):
    as_regex = re.compile("[IL]")
    psm_anno_dic = {}
    psm_anno_dic["anno"] = []
    psm_anno_dic["not_anno"] = []
    for scan, psms in scan_pep_map.items():
        for psm in psms:
            protein = psm["proteins"]
            if ("DECOY_" in protein) or (len(protein.split(",")) > 1):
                continue
            if orf_cds_dic[protein][3] == "intergenic":
                continue
            if orf_cds_dic[protein][1] in ["splice", "error_anno"]:
                continue
            pep_seq = psm["pep"]
            protein_seq = frame6[protein].seq
            org_start = int(protein.split("|")[4].split("-")[0])
            org_stop = int(protein.split("|")[4].split("-")[1])
            strand = protein.split("|")[3]
            contig = protein.split('|')[1]

            pep_seq_regex = as_regex.sub("[IL]", pep_seq)
            pep_starts = [m.start() for m in re.finditer(pep_seq_regex,
                                                         str(protein_seq))]
            if(pep_starts == []):
                print("Peptide not in protein")
                print(psm)
                print(pep_seq)
                raise
            if len(pep_starts) != 1:
                print("Multi hit in same protein!")
                print(psm)
                print(pep_seq)
                print(protein_seq)

            pep_start = pep_starts[0]
            if strand == "-1":
                pep_stop = org_stop - pep_start * 3
                pep_start = pep_stop - len(pep_seq) * 3
            else:
                pep_start = org_start + pep_start * 3
                pep_stop = pep_start + len(pep_seq) * 3

            genes = IntervalTree_dic[contig][pep_start:pep_stop]

            if len(genes) != 0:
                # print("Overlap")
                # print(psm)
                # print("{}-{}".format(pep_start, pep_stop))
                if orf_cds_dic[protein][0] != "None":
                    psm_anno_dic["anno"].append(psm[value])
                else:
                    psm_anno_dic["not_anno"].append(psm[value])
    psm_anno_dic["not_anno"].sort()
    psm_anno_dic["anno"].sort()
    return psm_anno_dic


def collect_genome(scan_pep_map, orf_cds_dic, frame_dic, value="e-value"):
    as_regex = re.compile("[IL]")
    psm_genome_dic = {}
    psm_genome_dic["anno"] = []
    psm_genome_dic["not_anno"] = []
    for scan, psms in scan_pep_map.items():
        for psm in psms:
            orf = psm["proteins"]
            if ("DECOY_" in orf) or (len(orf.split(",")) > 1):
                continue
            if orf_cds_dic[orf][0] == "None":
                psm_genome_dic["not_anno"].append(psm[value])
                continue
            try:
                start_anno = int(orf_cds_dic[orf][1]) / 3
            except ValueError:
                continue
            protein_seq = frame_dic[orf].seq
            pep_seq = psm["pep"]
            pep_seq_regex = as_regex.sub("[IL]", pep_seq)
            pep_stops = [m.start() + len(pep_seq) for m in
                         re.finditer(pep_seq_regex, str(protein_seq))]
            if any([pep_stop > start_anno for pep_stop in pep_stops]):
                psm_genome_dic["anno"].append(psm[value])
            else:
                psm_genome_dic["not_anno"].append(psm[value])
    psm_genome_dic["not_anno"].sort()
    psm_genome_dic["anno"].sort()
    return psm_genome_dic


def merge_results(res_dic, exp_dic):
    table = []
    for exp, values in res_dic.items():
        print(exp)
        scan_pep_map = exp_dic[exp]
        all_e_val = list(set(values["genome"]["not_anno"] +
                             values["genome"]["anno"] +
                             values["anno"]["not_anno"] +
                             values["anno"]["anno"]))
        all_e_val.sort()
        FDR_genome_interval = 0
        for e_val in all_e_val:
            # genome
            num_no_anno_above = bisect_left(values["genome"]["not_anno"],
                                            e_val)
            num_anno_above = bisect_left(values["genome"]["anno"], e_val)
            FDR_genome = round(float(num_no_anno_above) / float(num_anno_above +
                                                                num_no_anno_above
                                                                + 1) * 100.0, 2)

            if FDR_genome < FDR_genome_interval:
                continue
            FDR_genome_interval += 1

            for scan, psms in scan_pep_map.items():
                for psm in psms:
                    if psm["e-value"] == e_val:
                        FDR_decoy = psm["FDR"]

            num_no_anno_above = bisect_left(values["anno"]["not_anno"],
                                            e_val)
            num_anno_above = bisect_left(values["anno"]["anno"], e_val)

            FDR_anno = round((float(num_no_anno_above) / float(num_anno_above +
                                                               num_no_anno_above
                                                               + 1) * 1.2) *
                             100.0, 2)

            table.append("\t".join(map(str, [e_val, FDR_genome, FDR_anno,
                                             FDR_decoy, exp])))

    return table


def get_fdr_genomic(res_dic, FDR_decoy=1, category="anno"):
    num_anno_above = 0
    num_no_anno_above = 0
    for exp, values in res_dic.items():
        print(exp)
        num_no_anno_above += bisect_left(values[category]["not_anno"], FDR_decoy)
        num_anno_above += bisect_left(values[category]["anno"], FDR_decoy)

    if category == "anno":
        FDR = round((float(num_no_anno_above) / float(num_anno_above +
                                                      num_no_anno_above + 1) *
                     1.2) * 100.0, 2)
    else:
        FDR = round(float(num_no_anno_above) / float(num_anno_above +
                                                     num_no_anno_above + 1) *
                    100.0, 2)
    print(num_no_anno_above)
    print(num_anno_above)
    return FDR


def main():
    frame6 = SIHUMI.new_SIHUMI_data.load_6frame()
    orf_cds_dic = SIHUMI.new_SIHUMI_data.load_orf_cds_map()
    bed_file = SIHUMI.genome_dir.new_SIHUMI + "/ecoli.bed"
    crap_headers = SIHUMI.get_crap_headers()
    data_dir = "data/ecoli/"
    psm_base_dir = "/scr/k61san2/john/MS-Spektren/psm/6frame/ecoli/"
    e_cut_off = 10
    res_file = data_dir + "fdr_comp_new.tsv"

    IntervalTree_dic = {}
    with open(bed_file, "r") as file_handler:
        IntervalTree_dic = bed2IntervalTree_dic(file_handler, IntervalTree_dic)

    # collect ms
    exp_dic = {}
    '''
    exp_dic = collect_exp(psm_base_dir, crap_headers, e_cut_off)
    with open(data_dir + "big_exp_dic", "w") as file_handler:
        json.dump(exp_dic, file_handler)
    '''

    with open(data_dir + "big_exp_dic_new.json", "r") as file_handler:
        exp_dic = json.load(file_handler)

    '''
    result_dic = {}
    for experiment, scan_pep_map in exp_dic.items():
        result_dic[experiment] = {}
        result_dic[experiment]["anno"] = collect_anno(scan_pep_map,
                                                      IntervalTree_dic,
                                                      orf_cds_dic, frame6)

        result_dic[experiment]["genome"] = collect_genome(scan_pep_map,
                                                          orf_cds_dic, frame6)

    with open(data_dir + "fdr_result_dic.json", "w") as file_handler:
        json.dump(result_dic, file_handler)

    with open(data_dir + "fdr_result_dic.json", "r") as file_handler:
        result_dic = json.load(file_handler)

    with open(res_file, "w") as file_handler:
        file_handler.write("\n".join(merge_results(result_dic, exp_dic)))
    '''

    result_fdr_dic = {}
    for experiment, scan_pep_map in exp_dic.items():
        result_fdr_dic[experiment] = {}
        result_fdr_dic[experiment]["anno"] = collect_anno(scan_pep_map,
                                                          IntervalTree_dic,
                                                          orf_cds_dic, frame6,
                                                          value="FDR")

        result_fdr_dic[experiment]["genome"] = collect_genome(scan_pep_map,
                                                              orf_cds_dic,
                                                              frame6,
                                                              value="FDR")

    with open(data_dir + "result_dic_over_FDR.json", "w") as file_handler:
        json.dump(result_fdr_dic, file_handler)

    with open(data_dir + "result_dic_over_FDR.json", "r") as file_handler:
        result_fdr_dic = json.load(file_handler)


if __name__ == "__main__":
    main()

'''
import importlib
import fdr_anno

import json
import SIHUMI
from intervaltree import IntervalTree
import re
import fdr_anno
from bisect import bisect_left

frame6 = SIHUMI.new_SIHUMI_data.load_6frame()
orf_cds_dic = SIHUMI.new_SIHUMI_data.load_orf_cds_map()
bed_file = SIHUMI.genome_dir.new_SIHUMI + "/ecoli.bed"
crap_headers = SIHUMI.get_crap_headers()
data_dir = "data/ecoli/"
psm_base_dir = "/scr/k61san2/john/MS-Spektren/psm/6frame/ecoli/"
e_cut_off = 10
res_file = data_dir + "/fdr_comp.tsv"

IntervalTree_dic = {}
with open(bed_file, "r") as file_handler:
    IntervalTree_dic = fdr_anno.bed2IntervalTree_dic(file_handler, IntervalTree_dic)

# collect ms
exp_dic = {}

with open(data_dir + "big_exp_dic_new.json", "r") as file_handler:
    exp_dic = json.load(file_handler)

result_dic = {}

with open(data_dir + "fdr_result_dic.json", "r") as file_handler:
    result_dic = json.load(file_handler)

fdr_anno.merge_results(result_dic, exp_dic)
fdr_anno = importlib.reload(fdr_anno); print("\n".join(fdr_anno.merge_results(result_dic, exp_dic)))

values = result_dic["Anne_Mendler_P-1"]
e_val = 1
num_no_anno_above = bisect_left(values["anno"]["not_anno"], e_val)
num_anno_above = bisect_left(values["anno"]["anno"], e_val)

FDR_anno = round((float(num_no_anno_above) / float(num_anno_above + num_no_anno_above + 1) * 1.2) * 100.0, 2)

result_dic = {}
for experiment, scan_pep_map in exp_dic.items():
    result_dic[experiment] = {}
    result_dic[experiment]["anno"] = fdr_anno.collect_anno(scan_pep_map,
                                                           IntervalTree_dic,
                                                           orf_cds_dic, frame6)
    result_dic[experiment]["genome"] = fdr_anno.collect_genome(scan_pep_map,
                                                               orf_cds_dic,
                                                               frame6)

with open(data_dir + "result_dic_over_FDR.json", "r") as file_handler:
    result_dic = json.load(file_handler)
'''
