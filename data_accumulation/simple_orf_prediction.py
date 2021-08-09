#!/usr/bin/python3

import json
import re
from Bio.Data import CodonTable
from Bio import SeqIO
import os


def orf_prediction(orf_to_CDS, prot_dic, proteom_dic, genome, frame_dic):
    translation_table = CodonTable.ambiguous_dna_by_id[11]

    start_regex = re.compile('|'.join(translation_table.start_codons))
    canno_start_regex = re.compile("ATG|ATB|ATD|ATH|ATK|ATM|ATN|ATR|ATS|ATV|ATX|" +
                                   "DTG|MTG|NTG|RTG|VTG|WTG|XTG")
    as_regex = re.compile("[IL]")

    list_start_diff = []
    earlier_starts = {}
    for protein, psms in prot_dic["ecoli"]["6frame"].items():
        species = protein.split("|")[0]
        if str(orf_to_CDS[species][protein][1]) in ["None", "annotation error"]:
            continue

        len_protein_proteom = len(proteom_dic[orf_to_CDS[species][protein][0]])
        start_anno = orf_to_CDS[species][protein][1] / 3
        strand = protein.split("|")[3]
        start_protein = int(protein.split("|")[4].split("-")[0])
        stop_protein = int(protein.split("|")[4].split("-")[1])
        protein_seq = frame_dic[protein].seq
        chrom = protein.split("|")[1]

        if strand == "1":
            nuc_seq = genome[species][chrom][start_protein:stop_protein].seq
        else:
            nuc_seq = genome[species][chrom][start_protein:stop_protein].seq.reverse_complement()
        trans = nuc_seq.translate()
        if str(trans) != str(protein_seq):
            print("error! 6frame not equal to genome")
            print(protein)
            raise

        if "*" in trans:
            print("error! stop codon in translation")
            print(protein)
            raise
        pep_starts = []

        pep_start_dic = {}
        for psm in psms:
            pep_seq_regex = as_regex.sub("[IL]", psm["pep"])
            starts = [m.start() for m in re.finditer(pep_seq_regex,
                                                     str(protein_seq))]
            pep_start_dic[min(starts)] = psm["pep"]
            pep_starts += starts

        # start of startcodons
        canno_start_starts = []
        for m in re.finditer(canno_start_regex, str(nuc_seq)):
            if m.start() % 3 == 0:
                canno_start_starts.append(m.start())

        start_starts = []
        for m in re.finditer(start_regex, str(nuc_seq)):
            if m.start() % 3 == 0:
                start_starts.append(m.start())

        # calculate first start codon before mapped pep start
        min_start = min(pep_starts)
        # looks if a cannonical start exist if not take also alternatives,
        # if that does not exist 0 is assumed
        next_start_codon = max([start for start in canno_start_starts if start <=
                                min_start * 3] + [-1]) / 3
        if next_start_codon < 0:
            next_start_codon = max([start for start in start_starts if start <=
                                    min_start * 3] + [0]) / 3
        start_diff = next_start_codon - start_anno
        start_diff_ratio = start_diff / float(len_protein_proteom)
        '''
        if start_diff_ratio < 0:
            print(protein)
            print(protein_seq)
            print(start_anno)
            print(next_start_codon)
            print(start_diff)
            print(start_diff_ratio)
            print(len_protein_proteom)
            exit()
        '''

        if start_diff < 0:
            earlier_starts[protein] = psms
        list_start_diff.append([start_diff_ratio, len(psms), protein])

    return list_start_diff, earlier_starts


def load_data():
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    db_dir = data_dir + "/dbs"
    with open(db_dir + "/orf_to_CDS.json", "r") as file_handle:
        orf_to_CDS = json.load(file_handle)

    frame_dic = SeqIO.index(db_dir + "/SIHUMI_6frame.fasta", "fasta")

    proteom_dic = SeqIO.index(db_dir + "/SIHUMI_proteom.fasta", "fasta")

    genome_dir = data_dir + "/genome"
    genome = {}
    for genbank_file in os.listdir(genome_dir):
        if genbank_file.endswith(".gbk"):
            species = genbank_file.split(".")[0]
            genome_path = genome_dir + "/" + genbank_file
            genome[species] = SeqIO.to_dict(SeqIO.parse(genome_path, "genbank"))

    accu_data_dir = data_dir + "/accumulated_data"

    with open(accu_data_dir + "/prot_dic.json", "r") as file_handle:
        prot_dic = json.load(file_handle)
    return orf_to_CDS, prot_dic, proteom_dic, genome, frame_dic, accu_data_dir


def main():
    (orf_to_CDS, prot_dic, proteom_dic,
     genome, frame_dic, accu_data_dir) = load_data()

    list_start_diff, earlier_starts = orf_prediction(orf_to_CDS, prot_dic, proteom_dic, genome, frame_dic)

    with open(accu_data_dir + "/start_diffs", "w") as file_handle:
        print_start_diff = [str(entry[0]) + "," + str(entry[1]) for entry in list_start_diff]
        file_handle.write("\n".join(print_start_diff))

    with open(accu_data_dir + "/earlier_start.json", "w") as file_handle:
        json.dump(earlier_starts, file_handle)


if __name__ == "__main__":
    main()


'''
(orf_to_CDS, prot_dic, proteom_dic,
 genome, frame_dic, accu_data_dir) = load_data()
list_start_diff, earlier_starts = orf_prediction(orf_to_CDS, prot_dic, proteom_dic, genome, frame_dic)
with open(accu_data_dir + "/start_diffs", "w") as file_handle:
    print_start_diff = [str(entry[0]) + "," + str(entry[1]) for entry in list_start_diff]
    file_handle.write("\n".join(print_start_diff))

with open(accu_data_dir + "/earlier_start.json", "w") as file_handle:
    json.dump(earlier_starts, file_handle)
'''
