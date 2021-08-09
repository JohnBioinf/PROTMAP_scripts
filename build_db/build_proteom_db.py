#!/usr/bin/python3

import json
import os
from Bio import SeqIO


def make_proteom(data_dir):
    genome_dir = data_dir + "/genome/"
    SIHUMI_proteom_file = data_dir + "/dbs/SIHUMI_proteom.fasta"
    ecoli_proteom_file = data_dir + "/dbs/ecoli_proteom.fasta"

    SIHUMI_file_handle = open(SIHUMI_proteom_file, "w")
    ecoli_file_handle = open(ecoli_proteom_file, "w")
    for genbank_file in os.listdir(genome_dir):
        if genbank_file.endswith(".gbk"):
            genome = SeqIO.to_dict(SeqIO.parse(genome_dir + "/" +
                                               genbank_file, "genbank"))
            species = genbank_file.split(".")[0]
            for contig_name, contig in genome.items():
                for feature in contig.features:
                    if feature.type != "CDS":
                        continue
                    if "pseudo" in feature.qualifiers.keys():
                        continue
                    seq = str(feature.extract(genome[contig_name].seq).translate()[:-1])
                    seq.replace("*", "X")
                    if "protein_id" not in feature.qualifiers.keys():
                        protein_id = feature.qualifiers["locus_tag"][0]
                    else:
                        protein_id = feature.qualifiers["protein_id"][0]
                    strand = str(feature.location.strand)
                    start = int(feature.location.start)
                    stop = int(feature.location.end)
                    if strand == "-1":
                        start = start + 3
                    else:
                        stop = stop - 3
                    protein_id = (species + "|" + contig_name + "|" +
                                  protein_id + "|" + strand + "|" + str(start)
                                  + "-" + str(stop))
                    header = ">" + protein_id
                    SIHUMI_file_handle.write(header + "\n")
                    SIHUMI_file_handle.write(seq + "\n")
                    if species == "ecoli":
                        ecoli_file_handle.write(header + "\n")
                        ecoli_file_handle.write(seq + "\n")

    SIHUMI_file_handle.close()
    ecoli_file_handle.close()


def main():
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)["data_dir"]
    make_proteom(data_dir)
    return


if __name__ == "__main__":
    main()
