#!/usr/bin/python3

from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)

parser.add_argument("species", help="mini species name. E.g ecoli", type=str)
parser.add_argument("genome_dir", help="directory where the .gbk file resides",
                    type=str)
parser.add_argument("output_file", help="output bed file", type=str)

args = parser.parse_args()

print("loading genome")
genome = SeqIO.to_dict(SeqIO.parse(args.genome_dir + "/" + args.species + ".gbk", "genbank"))
print("finished loading genome")

rgb_colors_protein_dic = {}
# frame 1 red
rgb_colors_protein_dic[1] = "255,26,26"
# frame 2 blue
rgb_colors_protein_dic[2] = "26,26,255"
# frame 3 green
rgb_colors_protein_dic[3] = "45,185,45"
# frame 4 purple
rgb_colors_protein_dic[4] = "172,0,230"
# frame 5 orange
rgb_colors_protein_dic[5] = "230,115,0"
# frame 6 turquoise
rgb_colors_protein_dic[6] = "0,230,184"

rgb_colors_rna_dic = {}
# frame 1 red
rgb_colors_rna_dic[1] = "255,128,128"
# frame 2 blue
rgb_colors_rna_dic[2] = "128,128,255"
# frame 3 green
rgb_colors_rna_dic[3] = "153,230,153"
# frame 4 purple
rgb_colors_rna_dic[4] = "223,128,255"
# frame 5 orange
rgb_colors_rna_dic[5] = "255,191,128"
# frame 6 turquoise
rgb_colors_rna_dic[6] = "128,255,229"

# values for all
score = 0
track_itemRgb = "on"
blockCount = 1
blockStarts = 0
expCount = 1
expIds = 1
bed_list = []

for contig_name, contig in genome.items():
    track_name = track_description = contig_name
    track = "track name=\"{0}\" description=\"{1}\" itemRgb=\"{2}\"".format(track_name, track_description, track_itemRgb)
    bed_list.append(track)
    for feature in contig.features:
        if feature.type == "source":
            continue
        chromStart = thickStart = int(feature.location.start)
        chromEnd = thickEnd = int(feature.location.end)
        if(feature.location.strand > 0):
            strand = "+"
            frame = chromStart % 3 + 1
        else:
            frame = chromEnd % 3 + 4
            strand = "-"
        if feature.type == "CDS":
            if "protein_id" not in feature.qualifiers.keys():
                protein_id = feature.qualifiers["locus_tag"][0]
            else:
                protein_id = feature.qualifiers["protein_id"][0]
            if "pseudo" in feature.qualifiers.keys() or "join" in str(feature.location):
                continue
            name = "anno|" + contig_name + "|" + feature.type + "|" + protein_id.replace(" ", "").replace("|", ";")
            name = (name[:100] + '..') if len(name) > 100 else name
            itemRgb = rgb_colors_protein_dic[frame]
        elif "RNA" in feature.type:
            if "product" in feature.qualifiers.keys():
                product = feature.qualifiers["product"][0].replace("|", " ").replace(" ", "_")
            elif "gene" in feature.qualifiers.keys():
                product = feature.qualifiers["gene"][0].replace("|", " ").replace(" ", "_")
            elif "locus_tag" in feature.qualifiers.keys():
                product = feature.qualifiers["locus_tag"][0].replace("|", " ").replace(" ", "_")
            else:
                product = feature.type.replace("|", " ").replace(" ", "_")
            name = "anno|" + contig_name + "|" + feature.type + "|" + product
            name = (name[:100] + '..') if len(name) > 100 else name
            itemRgb = rgb_colors_rna_dic[frame]
        elif feature.type == "gene":
            continue
        elif feature.type == "assembly_gap":
            continue
        elif feature.type == "repeat_region":
            continue
        elif feature.type == "mobile_element":
            continue
        elif feature.type == "regulatory":
            continue
        elif feature.type == "rep_origin":
            continue
        elif "misc" in feature.type:
            continue
        elif feature.type == "gap":
            continue
        else:
            print("error unknown type")
            print(feature)
            exit(1)
        blockSizes = chromEnd - chromStart
        line = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13}".format(contig_name, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts, expCount, expCount)
        bed_list.append(line)

with open(args.output_file, "w") as the_file:
    the_file.write("\n".join(bed_list) + "\n")
