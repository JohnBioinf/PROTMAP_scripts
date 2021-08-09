#!/usr/bin/python3

import json
from Bio import SeqIO
import os
from glob import glob
from intervaltree import IntervalTree
import difflib
from math import log
import sys

'''
num_good_psm

num_exp

num_psm

species

contig

strand

start_orf

stop_orf

seq_orf

seq_len_orf

start

stop

seq

seq_len

start_genomic_context

stop_genomic_context

ucsc_link

name

num

candidate

blast_res

score_trans

genomic_contex_category

genomic_contex_strand

genomic_expression

'''


def mean_top_psms(psms):
    psm_scores = [psm["e-value"] for psm in psms]
    psm_scores.sort()
    return sum(psm_scores[0:3])/3


def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
    return s1[pos_a:pos_a+size]


def count_unique_psm(psms):
    seqs = list(set([psm["pep"] for psm in psms]))
    i = 0
    while i < len(seqs):
        seq = seqs[i]
        j = i + 1
        while j < len(seqs):
            overlap = get_overlap(seq, seqs[j])
            full_length = len(seq) + len(seqs[j]) - len(overlap)
            percent = len(overlap) / full_length * 100.0
            # print(seq)
            # print(seqs[j])
            # print(percent)
            # print()
            if percent > 80:
                del(seqs[j])
            else:
                j += 1
        i += 1
    # print("\n".join(seqs))
    return len(seqs)


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


def make_intervaltree(bed_dir):
    IntervalTree_dic = {}
    for bed_file in glob(bed_dir + "/*.bed"):
        with open(bed_file, "r") as file_handler:
            IntervalTree_dic = bed2IntervalTree_dic(file_handler, IntervalTree_dic)

    return IntervalTree_dic


def collect_general_results(info, candidate, psms, cand_bed_dir, genome_dir,
                            frame_dic, genome, i, selection_file, url_temp):
    info["name"] = selection_file + "_" + str(i)
    info["rank"] = i
    info["candidate"] = candidate

    info["num_good_psm"] = len([1 for psm in psms
                                if psm["num"] == "1" and
                                len(psm["proteins"].split(",")) == 1])
    info["num_exp"] = len(set([psm["experiment"] for psm in psms]))
    info["num_psm"] = len(psms)
    info["num_good_psm"] = len([1 for psm in psms if psm["num"] == "1" and
                                len(psm["proteins"].split(","))])
    info["mean_top_psms"] = log(mean_top_psms(psms), 10) * -1
    info["unique_psms"] = count_unique_psm(psms)
    info["species"] = candidate.split("|")[0]
    info["contig"] = candidate.split("|")[1]
    info["strand"] = candidate.split("|")[3]

    info["start_orf"] = candidate.split("|")[4].split("-")[0]
    info["stop_orf"] = candidate.split("|")[4].split("-")[1]
    info["seq_orf"] = str(frame_dic[candidate].seq)
    info["seq_len_orf"] = len(info["seq_orf"])

    with open(cand_bed_dir.format(info["species"]), "r") as file_handler:
        for line in file_handler:
            if info["name"] + " " in line:
                info["start"] = line.split(" ")[1]
                info["stop"] = line.split(" ")[2]

    nuc_seq = genome[info["species"]][info["contig"]][int(info["start"]):int(info["stop"])].seq
    if info["strand"] == "1":
        info["seq"] = str(nuc_seq.translate())
        info["start_codon"] = str(nuc_seq[0:3])
    else:
        info["start_codon"] = str(nuc_seq[0:3])
        info["seq"] = str(nuc_seq.reverse_complement().translate())
    info["seq_len"] = len(info["seq"])

    len_contig = len(genome[info["species"]][info["contig"]].seq)
    start_diffs = [0]
    stop_diffs = []
    for feature in genome[info["species"]][info["contig"]].features:
        if feature.type != "CDS":
            continue
        if "pseudo" in feature.qualifiers.keys():
            continue
        if "join" in str(feature.location):
            continue
        gene_start = feature.location.start
        gene_stop = feature.location.end
        start_diffs.append(int(info["start"]) - int(gene_start))
        stop_diffs.append(int(gene_stop) - int(info["stop"]))

    info["start_genomic_context"] = (int(info["start"]) -
                                     min([s for s in start_diffs if s > 0]
                                         + [len_contig]))
    info["stop_genomic_context"] = (min([s for s in stop_diffs if s > 0]
                                        + [0]) + int(info["stop"]))

    info["ucsc_link"] = url_temp.format(info["species"], info["contig"],
                                        info["start"], info["stop"])
    return info


def check_hypo_blast_res(name):
    key_names = ["ypothetic", "utativ", "ncharacterised", "predicted",
                 "ncharacterized", "conserved domain protein", "redicted",
                 "unknown"]
    for key_name in key_names:
        if key_name in name:
            return True
    return False


def collect_blast_res(info, blast_file, e_val_cutoff=float("1.0e-5")):
    blast_res = []
    blast_category = ""
    with open(blast_file + ".fasta.tsv", "r") as file_handler:
        for line in file_handler:
            line = line.split("\t")
            e_value = float(line[1])
            if e_value > e_val_cutoff:
                break
            titels = line[5].split("<>")
            if blast_category == "":
                if not any([check_hypo_blast_res(title) for title in titels]):
                    for title in titels:
                        if not check_hypo_blast_res(title):
                            blast_category = title
                            break
            query_cov = line[7]
            ref_seq = line[4]
            seq_sim = line[3]
            blast_res.append([titels, e_value, query_cov, seq_sim, ref_seq])
    if blast_res == []:
        blast_category = "Novel"
    elif blast_category == "":
        blast_category = "Hypothetical"

    info["blast_res"] = blast_res
    info["blast_category"] = blast_category
    return info


def quantile(l, q):
    l.sort()
    return l[int(float(len(l)) * q)]


def load_transcriptom_data(coverage_dic_file, background_dic_file, trans_file):
    if os.path.exists(coverage_dic_file):
        with open(coverage_dic_file, "r") as file_handler:
            coverage_dic = (json.load(file_handler))

        with open(background_dic_file, "r") as file_handler:
            background_dic = (json.load(file_handler))

    else:
        coverage_dic = {}
        with open(trans_file, "r", newline="") as coverage_file:
            for line in coverage_file:
                line = line.replace("\n", "").split(" ")
                contig = line[0]
                start = int(line[1])
                cov = float(line[3])
                if contig in coverage_dic:
                    coverage_dic[contig].append([start, cov])
                else:
                    coverage_dic[contig] = [[start, cov]]

        background_dic = {}
        for contig, coverage in coverage_dic.items():
            background_dic[contig] = quantile([cov[1] for cov in coverage], 0.5)

        with open(background_dic_file, "w") as file_handler:
            file_handler.write(json.dumps(background_dic))

        with open(coverage_dic_file, "w") as file_handler:
            file_handler.write(json.dumps(coverage_dic))

    return coverage_dic, background_dic


def collect_transcriptom(info, candidate, coverage_dic, background_dic):
    contig = candidate.split("|")[1]
    start = int(candidate.split("|")[4].split("-")[0])
    stop = int(candidate.split("|")[4].split("-")[1])
    background_noise = background_dic[contig]
    cov_list = []
    for line in coverage_dic[contig]:
        if line[0] >= start:
            if line[0] > stop:
                break
            cov_list.append(line)
    if len(cov_list) == 0:
        score_trans = 0
    else:
        over_background_list = [0 if cov[1] < background_noise else 1 for cov in cov_list]
        score_trans = round(sum(over_background_list)/len(over_background_list), 2)
    info["score_trans"] = score_trans
    return info


def collect_genomic_context(info, candidate, bed_intervaltree, CDS_to_orf,
                            prot_dic):
    contig = candidate.split("|")[1]
    strand = candidate.split("|")[3]
    start = int(candidate.split("|")[4].split("-")[0])
    stop = int(candidate.split("|")[4].split("-")[1])
    gene_intervals = bed_intervaltree[contig][start:stop]
    translation = 0
    if len(gene_intervals) == 0:
        info["genomic_context_category"] = "no_overlap"
        info["genomic_context_strand"] = "None"
    else:
        info["genomic_context_category"] = "overlap"
        info["genomic_context_strand"] = "other_strand"

    for interval in gene_intervals:
        gene = interval.data.split("|")[-1]
        for CDS in CDS_to_orf[info["species"]].keys():
            if gene in CDS:
                orfs = CDS_to_orf[info["species"]][CDS][0]

        for orf in orfs:
            gene_strand = orf.split("|")[3]
            if gene_strand == strand:
                info["genomic_context_strand"] = "same_strand"
            if orf not in prot_dic["SIHUMI"]["6frame"]:
                continue
            psms = prot_dic["SIHUMI"]["6frame"][orf]

            psm_context = [psm for psm in psms if psm["num"] == 1 and
                           len(psm["proteins"].split(",")) == 1 and
                           psm["FDR"] < 1]
            num_unique_psms = count_unique_psm(psm_context)
            translation += num_unique_psms
    info["genomic_context_translated"] = translation
    return info


def classification(info):
    decision_process = "{num(PSMs) > 5} ->"
    if info["mean_top_psms"] > 3.5:
        decision_process += "{s^ > 3.5} ->"
        if(info["genomic_context_strand"] == "same_strand" and
           info["genomic_context_translated"] > 6):
            decision_process += "{translation on same strand (frameshift candidate)}"
            classification = "True candidate"
        else:
            decision_process += "{no translation on same strand} ->"
            if info["unique_psms"] > 2:
                decision_process += "{num(unique peptides)} > 2"
                classification = "True candidate"
            else:
                decision_process += "{num(unique peptides)} < 2"
                if info["score_trans"] > 0.7:
                    decision_process += "{transkription is high}"
                    classification = "True candidate"
                else:
                    decision_process += "{transkription is low}"
                    classification = "False candidate"
    else:
        decision_process += "{s^ < 3.5} ->"
        if info["unique_psms"] > 3:
            decision_process += "{num(unique peptides) > 3}"
            classification = "True candidate"
        else:
            decision_process += "{num(unique peptides) < 3}"
            if info["mean_top_psms"] > 2.5:
                decision_process += "{s^ > 2.5} ->"
                if info["genomic_context_category"] == "no_overlap":
                    decision_process += "{no overlap with annotation}"
                    classification = "True candidate"
                else:
                    decision_process += "{overlap with annotation}"
                    classification = "False candidate"
            else:
                decision_process += "{s^ < 2.5}"
                classification = "False candidate"

    info["classification"] = classification
    info["decision_process"] = decision_process
    return info


def make_info(cand_bed_dir, genome_dir, frame_dic, genome, selection_file,
              url_temp, cand_list, prot_dic, coverage_dic, background_dic,
              CDS_to_orf, bed_intervaltree, blast_dir):
    info_dic = {}

    for i, candidate in enumerate(cand_list):
        psms = prot_dic["SIHUMI"]["6frame"][candidate]
        info = {}
        info = collect_general_results(info, candidate, psms, cand_bed_dir,
                                       genome_dir, frame_dic, genome, i + 1,
                                       selection_file, url_temp)
        info = collect_blast_res(info, blast_dir + info["name"])
        info = collect_transcriptom(info, candidate, coverage_dic,
                                    background_dic)
        info = collect_genomic_context(info, candidate, bed_intervaltree,
                                       CDS_to_orf, prot_dic)
        info = classification(info)
        info_dic[candidate] = info

    return info_dic


def main():
    selection_cut_off = int(sys.argv[1])
    selection_file = "nov_psm" + str(selection_cut_off)

    with open("./parameters.json", "r") as file_handle:
        json_content = json.load(file_handle)
        data_dir = json_content['data_dir']
        hub_id = json_content['hub_id']
        session_id = json_content['session_id']

    db_dir = data_dir + "/dbs"
    cand_dir = data_dir + "/candidates"
    genome_dir = data_dir + "/genome"
    trans_dir = data_dir + "/transcriptom"
    accu_data_dir = data_dir + "/accumulated_data"

    frame_dic = SeqIO.index(db_dir + "/SIHUMI_6frame.fasta", "fasta")

    genome = {}
    for genbank_file in os.listdir(genome_dir):
        if genbank_file.endswith(".gbk"):
            species = genbank_file.split(".")[0]
            genome_path = genome_dir + "/" + genbank_file
            genome[species] = SeqIO.to_dict(SeqIO.parse(genome_path, "genbank"))

    with open(db_dir + "/CDS_to_orf.json", "r") as file_handle:
        CDS_to_orf = json.load(file_handle)

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)

    bed_dir = cand_dir + "/bed_dir/" + selection_file + "/{}_cand.bed"
    blast_dir = cand_dir + "/blast/{}/".format(selection_file)

    url_temp_base = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hub_{0}_{2}&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position={2}%3A{2}-{2}{1}"
    url_temp = url_temp_base.format(hub_id, session_id, "{}")
    with open(cand_dir + "/" + selection_file + "_list.json", "r") as file_handler:
        cand_list = json.load(file_handler)

    coverage_dic_file = trans_dir + "/coverage_dic.json"
    background_dic_file = trans_dir + "/background_dic.json"
    trans_file = trans_dir + "/average.wig"
    coverage_dic, background_dic = load_transcriptom_data(coverage_dic_file,
                                                          background_dic_file,
                                                          trans_file)

    bed_intervaltree = make_intervaltree(genome_dir)
    info_dic = make_info(bed_dir, genome_dir, frame_dic, genome,
                         selection_file, url_temp, cand_list, prot_dic,
                         coverage_dic, background_dic, CDS_to_orf,
                         bed_intervaltree, blast_dir)

    with open(cand_dir + "/" + selection_file + "_info_dic.json", "w") as file_handler:
        json.dump(info_dic, file_handler)


if __name__ == "__main__":
    main()


'''
selection_cut_off = 6
selection_file = "nov_psm" + str(selection_cut_off)

with open("./parameters.json", "r") as file_handle:
    json_content = json.load(file_handle)
    data_dir = json_content['data_dir']
    hub_id = json_content['hub_id']
    session_id = json_content['session_id']

db_dir = data_dir + "/dbs"
cand_dir = data_dir + "/candidates"
genome_dir = data_dir + "/genome"
trans_dir = data_dir + "/transcriptom"
accu_data_dir = data_dir + "/accumulated_data"

frame_dic = SeqIO.index(db_dir + "/SIHUMI_6frame.fasta", "fasta")

genome = {}
for genbank_file in os.listdir(genome_dir):
    if genbank_file.endswith(".gbk"):
        species = genbank_file.split(".")[0]
        genome_path = genome_dir + "/" + genbank_file
        genome[species] = SeqIO.to_dict(SeqIO.parse(genome_path, "genbank"))

with open(db_dir + "/CDS_to_orf.json", "r") as file_handle:
    CDS_to_orf = json.load(file_handle)

with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
    prot_dic = json.load(file_handle)

bed_dir = cand_dir + "/bed_dir/" + selection_file + "/{}_cand.bed"
blast_dir = cand_dir + "/blast/{}/".format(selection_file)

url_temp_base = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hub_{0}_{2}&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position={2}%3A{2}-{2}{1}"
url_temp = url_temp_base.format(hub_id, session_id, "{}")
with open(cand_dir + "/" + selection_file + "_list.json", "r") as file_handler:
    cand_list = json.load(file_handler)

coverage_dic_file = trans_dir + "/coverage_dic.json"
background_dic_file = trans_dir + "/background_dic.json"
trans_file = trans_dir + "/average.wig"
coverage_dic, background_dic = load_transcriptom_data(coverage_dic_file,
                                                      background_dic_file,
                                                      trans_file)

bed_intervaltree = make_intervaltree(genome_dir)
info = {}

info_dic = make_info(bed_dir, genome_dir, frame_dic, genome, selection_file, url_temp, cand_list, prot_dic, coverage_dic, background_dic, CDS_to_orf, bed_intervaltree, blast_dir)

with open(cand_dir + "/" + selection_file + "_info_dic.json", "w") as file_handler:
    json.dump(info_dic, file_handler)

#
candidate = "bact|CP040530.1|512796|-1|310719-311007"
i = 24
psms = prot_dic["SIHUMI"]["6frame"][candidate]
info = collect_general_results(info, candidate, psms, bed_dir, genome_dir, frame_dic, genome, i + 1, selection_file, url_temp)
info = collect_genomic_context(info, candidate, bed_intervaltree, CDS_to_orf, prot_dic)

psm_context = [psm for psm in psms if psm["num"] == 1 and len(psm["proteins"].split(",")) == 1 and psm["FDR"] < 1]

'''
