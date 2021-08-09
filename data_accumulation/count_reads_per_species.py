#!/usr/bin/python3

import json


def load_data():
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    with open("./SIHUMI_info_dic.json", "r") as file_handle:
        SIHUMI_info_dic = json.load(file_handle)

    trans_dir = data_dir + "/transcriptom"
    genome_dir = data_dir + "/genome"
    accu_data_dir = data_dir + "/accumulated_data"

    return trans_dir, SIHUMI_info_dic, accu_data_dir, genome_dir


def count_reads_4_species(trans_dir):
    read_count_dic = {}
    with open(trans_dir + "/average.wig", "r") as file_handle:
        for line in file_handle:
            line = line[:-1].split(" ")
            contig = line[0]
            reads = int(line[3])
            if contig not in read_count_dic:
                read_count_dic[contig] = 0
            read_count_dic[contig] += reads

    return read_count_dic


def get_size_of_genome(genome_dir):
    genome_size_dic = {}
    with open(genome_dir + "/chrom.sze", "r") as file_handle:
        for line in file_handle:
            line = line[:-1].split("\t")
            genome_size_dic[line[0]] = int(line[-1])

    return genome_size_dic


def main():
    trans_dir, SIHUMI_info_dic, accu_data_dir, genome_dir = load_data()

    read_count_dic = count_reads_4_species(trans_dir)

    genome_size_dic = get_size_of_genome(genome_dir)

    # normalize
    res = ""
    for species, info_dic in SIHUMI_info_dic.items():
        species_norm = 0.0
        for contig in info_dic["contigs"]:
            contig_sum = [r for c, r in read_count_dic.items() if c == contig]
            if contig_sum == []:
                contig_sum = [0]
            species_norm += contig_sum[0] / genome_size_dic[contig]

        res += info_dic["short_name"] + "\t" + str(species_norm) + "\n"

    with open(accu_data_dir + "/reads_species_ratio.tsv", "w") as file_handler:
        file_handler.write(res)


if __name__ == "__main__":
    main()

'''
trans_dir, SIHUMI_info_dic, accu_data_dir, genome_dir = load_data()

read_count_dic = count_reads_4_species(trans_dir)
genome_size_dic = get_size_of_genome(genome_dir)

res = ""
for species, info_dic in SIHUMI_info_dic.items():
    species_norm = 0.0
    for contig in info_dic["contigs"]:
        contig_sum = [r for c, r in read_count_dic.items() if c == contig]
        if contig_sum == []:
            contig_sum = [0]
        species_norm += contig_sum[0] / genome_size_dic[contig]

    res += info_dic["short_name"] + "\t" + str(species_norm) + "\n"
'''
