#!/usr/bin/python3

import json
import re
from Bio.Data import CodonTable
from Bio import SeqIO
import sys
import os


def make_bed_and_fasta(cand_list, frame_dic, genome, prot_dic, selection_file):
    bed_dic = {}
    fasta_dic = {}

    frame_color_dic = {}

    frame_color_dic[1] = "255,26,26"
    frame_color_dic[2] = "26,26,255"
    frame_color_dic[3] = "45,185,45"
    frame_color_dic[4] = "172,0,230"
    frame_color_dic[5] = "230,115,0"
    frame_color_dic[6] = "0,230,184"

    bed_temp = "{0} {1} {2} {3} 0 {4} {1} {2} {5} 1 {6} 0 1 1"
    track_description = ('track name="Candidates" description="Candidate"'
                         'visibility=2 itemRgb="On"')
    translation_table = CodonTable.ambiguous_dna_by_id[11]

    start_regex = re.compile('|'.join(translation_table.start_codons))
    canno_start_regex = re.compile("ATG|ATB|ATD|ATH|ATK|ATM|ATN|ATR|ATS|ATV|" +
                                   "ATX|DTG|MTG|NTG|RTG|VTG|WTG|XTG")
    as_regex = re.compile("[IL]")
    for i, orf in enumerate(cand_list):
        species = orf.split("|")[0]
        psms = prot_dic["SIHUMI"]["6frame"][orf]

        strand = "+" if orf.split("|")[3] == "1" else "-"
        start = int(orf.split("|")[4].split("-")[0])
        stop = int(orf.split("|")[4].split("-")[1])
        chrom = orf.split("|")[1]
        protein_seq = frame_dic[orf].seq

        if strand == "+":
            nuc_seq = genome[species][chrom][start:stop].seq
        else:
            nuc_seq = genome[species][chrom][start:stop].seq.reverse_complement()
        trans = nuc_seq.translate()
        if str(trans) != str(protein_seq):
            print("error! 6frame not equal to genome")
            print(orf)
            raise

        if "*" in trans:
            print("error! stop codon in translation")
            print(orf)
            raise
        pep_starts = []

        pep_start_dic = {}
        for psm in psms:
            peptide = psm["pep"]
            pep_seq_regex = as_regex.sub("[IL]", peptide)
            starts = [m.start() for m in re.finditer(pep_seq_regex,
                                                     str(protein_seq))]
            pep_start_dic[min(starts)] = peptide
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
        '''
        for psm in psms:
            peptide = psm["pep"]
            if peptide == pep_start_dic[min_start]:
                print(peptide)
                break
        '''
        # looks if a cannonical start exist if not take also alternatives,
        # if that does not exist 0 is assumed
        next_start_codon = max([start for start in canno_start_starts if start <=
                                min_start * 3] + [-1])
        if next_start_codon < 0:
            next_start_codon = max([start for start in start_starts if start <=
                                    min_start * 3] + [0])
        '''
        print(next_start_codon)
        print(start_starts)
        print(canno_start_starts)
        print(min_start)
        '''
        if strand == "+":
            start = start + next_start_codon
            frame = start % 3 + 1
        else:
            stop = stop - next_start_codon
            frame = stop % 3 + 4
        rgb = frame_color_dic[frame]
        name = selection_file + "_" + str(i + 1)
        if strand == "+":
            cand_seq = genome[species][chrom][start:stop].seq.translate()
        else:
            cand_seq = genome[species][chrom][start:stop].seq.reverse_complement().translate()

        fasta_dic[name] = ">" + name + "\n" + str(cand_seq) + "\n"
        if species not in bed_dic.keys():
            bed_dic[species] = [track_description]

        bed_dic[species].append(bed_temp.format(chrom, start, stop, name,
                                                strand, rgb, stop - start))

    return bed_dic, fasta_dic


def main():
    selection_cut_off = int(sys.argv[1])
    selection_file = "nov_psm" + str(selection_cut_off)

    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    with open("./SIHUMI_info_dic.json", "r") as file_handle:
        SIHUMI_info_dic = json.load(file_handle)

    db_dir = data_dir + "/dbs"
    frame_path = db_dir + "/SIHUMI_6frame.fasta"
    frame_dic = SeqIO.index(frame_path, "fasta")

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

    cand_dir = data_dir + "/candidates"
    with open(cand_dir + "/" + selection_file + "_list.json", "r") as file_handler:
        cand_list = json.load(file_handler)

    bed_dic, fasta_dic = make_bed_and_fasta(cand_list, frame_dic, genome,
                                            prot_dic, selection_file)

    bed_dir = cand_dir + "/bed_dir/" + selection_file + "/"
    for species in SIHUMI_info_dic.keys():
        with open(bed_dir + species + "_cand.bed", "w") as file_handler:
            file_handler.write("\n".join(bed_dic[species]) + "\n")

    blast_dir = cand_dir + "/blast/{}/".format(selection_file)
    for name, content in fasta_dic.items():
        with open(blast_dir + name + ".fasta", "w") as file_handler:
            file_handler.write(content)


if __name__ == "__main__":
    main()


'''
'''
