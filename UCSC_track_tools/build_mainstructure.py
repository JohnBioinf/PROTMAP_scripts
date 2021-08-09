#!/usr/bin/python3

import os
import shutil
import subprocess
import json


def system_call(call_list):
    p = subprocess.Popen(call_list, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    p.wait()
    if p.returncode != 0:
        print("error:")
        print(p.communicate()[0].decode('utf-8'))
        exit()


def load_genome(genome_dir, out_dir, species):
    genome_commands = ["./UCSC_track_tools/load_genome.sh", genome_dir,
                       out_dir, species]
    system_call(genome_commands)


def load_annotation(output_dir, species, genome_dir):
    annotaion_commands = ["./UCSC_track_tools/load_annotation.sh", species, output_dir, tmp_dir,
                          genome_dir]
    system_call(annotaion_commands)


def load_massSpec(massSpec_dir, output_dir, species, db, tmp_dir, organism):
    massSpec_commands = ["./UCSC_track_tools/load_massSpec.sh", species, output_dir,
                         massSpec_dir, db, tmp_dir, organism]
    system_call(massSpec_commands)


def load_candidates(cand_dir, output_dir, species, tmp_dir):
    candidates_commands = ["./UCSC_track_tools/load_candidate.sh", cand_dir, species,
                           output_dir, tmp_dir]
    system_call(candidates_commands)


def clean_dir(directory):
    if(os.path.isdir(directory)):
        if(os.listdir(directory) != []):
            for f in os.listdir(directory):
                if os.path.isdir(directory + '/' + f):
                    shutil.rmtree(directory + '/' + f)
                else:
                    os.remove(directory + "/" + f)
    else:
        os.mkdir(directory)


# names files etc
with open("./parameters.json", "r") as file_handle:
    param = json.load(file_handle)
    data_dir = param['data_dir']
    tmp_dir = param['tmp_dir']

with open("./SIHUMI_info_dic.json", "r") as file_handle:
    SIHUMI_info_dic = json.load(file_handle)

genome_dir = data_dir + "/genome"
transcriptom_dir = data_dir + "/transcriptom"
track_dir = data_dir + "/SIHUMI_track_hub"
cand_dir = data_dir + "/candidates/bed_dir/nov_psm6"
ms_dir = data_dir + "/ms"

# combination of ms dirs that will be loaded. each combination will be on track
ms_combinations = [["6frame", "SIHUMI"], ["6frame", "ecoli"]]
# Transcriptom data
transcriptome_file = transcriptom_dir + "/average.bw"

# structure

# SIHUMI_track/
# +-- genomes.txt
# +-- groups.txt
# +-- hub.txt
# +-- anaero
# |   +-- anaero.2bit
# |   +-- anaero_annotation.html
# |   +-- anaero.chrom.sizes
# |   +-- anaero.jpg
# |   +-- bbi
# |   |   +-- anaero_annotation.bb
# |   |   +-- mapping_1_sorted.bigwig
# |   |   +-- mapping_2_sorted.bigwig
# |   |   +-- massSpec.bb
# |   |   +-- massSpec_SIHUMI_hannes.bb
# |   |   \-- massSpec_SIHUMI_steffi.bb
# |   +-- description.html
# |   +-- massSpec.html
# |   \-- trackDb.txt
# |
# ...


# templates
address_txt = ("http://www.bioinf.uni-leipzig.de/Publications/SUPPLEMENTS/"
               "20-002/SIHUMI_track_hub/hub.txt\n")
"http://www.bioinf.uni-leipzig.de/Publications/SUPPLEMENTS/20-002/SIHUMI_track_hub/hub.txt"
hub_txt = ("hub SIHUMI\n"
           "shortLabel SIHUMI\n"
           "longLabel simplified human intestinal microbiota HUB\n"
           "genomesFile genomes.txt\n"
           "descriptionUrl https://pubmed.ncbi.nlm.nih.gov/34039272/\n"
           "email john@bioinf.uni-leipzig.de\n")
groups_txt_template = "name {0}\nlabel {1}\npriority {2}\ndefaultIsClosed {3}"
genomes_txt = ("genome {0}\n"
               "trackDb {0}/trackDb.txt\n"
               "groups groups.txt\n"
               "description {1} \n"
               "twoBitPath {0}/{0}.2bit\n"
               "organism {1}\n"
               "defaultPos {4}:1-2000\n"
               "orderKey {3}\n"
               "scientificName {2}\n"
               "htmlPath {0}/description.html\n")

with open("UCSC_track_tools/description_html/massSpec_SIHUMI.html", "r") as file_handle:
    massSpec_SIHUMI_html = file_handle.read()

with open("UCSC_track_tools/description_html/massSpec_ecoli.html", "r") as file_handle:
    massSpec_ecoli_html = file_handle.read()

annotation_html_dic = {}
for s in SIHUMI_info_dic:
    file_name = "UCSC_track_tools/description_html/annotation_" + s + ".html"
    with open(file_name, "r") as file_handle:
        annotation_html_dic[s] = file_handle.read()

hub_html_dic = {}
for s in SIHUMI_info_dic:
    file_name = "UCSC_track_tools/description_html/" + s + "_hub.html"
    with open(file_name, "r") as file_handle:
        hub_html_dic[s] = file_handle.read()

with open("UCSC_track_tools/description_html/transcriptome.html", "r") as file_handle:
    transcriptome_html = file_handle.read()

with open("UCSC_track_tools/description_html/candidates.html", "r") as file_handle:
    candidate_html = file_handle.read()

# tracks
# 1 annotation
annotation_trackDb_txt = ("track annotation\n"
                          "longLabel Annotation\n"
                          "shortLabel Annotation\n"
                          "priority 1\n"
                          "visibility pack\n"
                          "itemRgb on\n"
                          "bigDataUrl bbi/annotation.bb\n"
                          "type bigBed 12 +\n"
                          "group annotation\n"
                          "html annotation\n"
                          "searchIndex name\n")

# 2 candidate
candidate_trackDb_txt = ("track candidate\n"
                         "longLabel Candidate\n"
                         "shortLabel Cand\n"
                         "priority 2\n"
                         "visibility pack\n"
                         "itemRgb on\n"
                         "bigDataUrl bbi/cand.bb\n"
                         "type bigBed 12 +\n"
                         "group candidate\n"
                         "html candidate\n")

# 3 transcriptome
transcriptome_trackDb_txt = ("track transcriptome\n"
                             "longLabel Transcriptome\n"
                             "shortLabel Trans\n"
                             "priority 3\n"
                             "visibility full\n"
                             "bigDataUrl bbi/average.bw\n"
                             "type bigWig\n"
                             "group transcriptome\n"
                             "html transcriptome\n")

# 4 mass spec
massSpec_trackDb_txt = ("track massSpec\n"
                        "longLabel Mass Spectrometry\n"
                        "shortLabel Mass Spec\n"
                        "priority 4\n"
                        "visibility squish\n"
                        "itemRgb on\n"
                        "bigDataUrl bbi/massSpec_6frame_SIHUMI.bb\n"
                        "type bigBed 12 +\n"
                        "group massSpec_SIHUMI\n"
                        "html massSpec_SIHUMI\n"
                        "searchIndex name\n")

# 5 mass spec only for ecoli
ecoli_massSpec_trackDb_txt = ("track massSpec_ecoli\n"
                              "longLabel Mass Spectrometry Ecoli\n"
                              "shortLabel Mass Spec Ecoli\n"
                              "priority 5\n"
                              "visibility squish\n"
                              "itemRgb on\n"
                              "bigDataUrl bbi/massSpec_6frame_ecoli.bb\n"
                              "type bigBed 12 +\n"
                              "group massSpec_ecoli\n"
                              "html massSpec_ecoli\n"
                              "searchIndex name\n")
trackDb_txt = "\n".join((annotation_trackDb_txt, candidate_trackDb_txt,
                         transcriptome_trackDb_txt, massSpec_trackDb_txt))

# infos for templates

# build tuple bundles to fit in templates with SIHUMI data
description_html_bundle = {}
genomes_txt_bundle = {}
for species in SIHUMI_info_dic.keys():
    genomes_txt_bundle[species] = (species,
                                   SIHUMI_info_dic[species]["short_name"],
                                   SIHUMI_info_dic[species]["long_name"], "{0}",
                                   SIHUMI_info_dic[species]["contigs"][0])

# groups
group_annotation_bundle = ("annotation", "Annotation", 1, 1)
group_cand_bundle = ("candidate", "Candidate", 2, 1)
group_transcriptome_bundle = ("transcriptome", "Transcriptome", 3, 1)
group_massSpec_bundle = ("massSpec", "Mass Spectrometry", 4, 1)
groups_bundle_list = [group_annotation_bundle, group_cand_bundle,
                      group_transcriptome_bundle, group_massSpec_bundle]

# ---------- Build Structure -------------

# build dir
print("Building clean file structure")
clean_dir(track_dir)
for s in SIHUMI_info_dic.keys():
    os.mkdir(track_dir + "/" + s)
    os.mkdir(track_dir + "/" + s + "/bbi")
print("Write basic files")

# print files
# hub.txt

with open(track_dir + "/hub.txt", "w") as file_handle:
    file_handle.write(hub_txt)

with open(track_dir + "/address.txt", "w") as file_handle:
    file_handle.write(address_txt)

# groups.txt
groups_txt = ""
for group_bundle in groups_bundle_list:
    groups_txt = groups_txt + groups_txt_template.format(*group_bundle) + "\n\n"
with open(track_dir + "/groups.txt", "w") as file_handle:
    file_handle.write(groups_txt[:-2])

# genomes.txt
new_genomes_txt = ""
i = 1
for s in SIHUMI_info_dic.keys():
    new_genomes_txt = new_genomes_txt + genomes_txt.format(*genomes_txt_bundle[s]) + "\n"
    new_genomes_txt = new_genomes_txt.format(i)
    i = i + 1
new_genomes_txt = new_genomes_txt[:-1]
with open(track_dir + "/genomes.txt", "w") as file_handle:
    file_handle.write(new_genomes_txt)

print("building species specific files")
for s in SIHUMI_info_dic.keys():
    print(s)
    with open(track_dir + "/" + s + "/description.html", "w") as file_handle:
        file_handle.write(hub_html_dic[s])

    if s == "ecoli":
        current_trackDb_txt = trackDb_txt + "\n" + ecoli_massSpec_trackDb_txt
    else:
        current_trackDb_txt = trackDb_txt

    print("load genome")
    load_genome(genome_dir, track_dir + "/" + s, s)
    print("load annotation")
    load_annotation(track_dir + "/" + s, s, genome_dir)
    print("load candidates")
    load_candidates(cand_dir, track_dir + "/" + s, s, tmp_dir)
    print("load transcriptome")
    shutil.copyfile(transcriptome_file, track_dir + "/" + s + "/bbi/" +
                    transcriptome_file.split("/")[-1])
    print("load ms")
    for ms_combination in ms_combinations:
        db = ms_combination[0]
        organism = ms_combination[1]
        if organism == "ecoli":
            if s != "ecoli":
                continue
            res_dir = (ms_dir + "/{}/bed/{}/").format(organism, db)
        else:
            res_dir = (ms_dir + "/{}/bed/{}/").format(organism, db)
        output_dir = track_dir + "/" + s
        load_massSpec(res_dir, output_dir, s, db, tmp_dir, organism)

    with open(track_dir + "/" + s + "/trackDb.txt", "w") as file_handle:
        file_handle.write(current_trackDb_txt.format(s))

    with open(track_dir + "/" + s + "/massSpec_SIHUMI.html", "w") as file_handle:
        file_handle.write(massSpec_SIHUMI_html)

    if s == "ecoli":
        with open(track_dir + "/" + s + "/massSpec_ecoli.html", "w") as file_handle:
            file_handle.write(massSpec_ecoli_html)

    with open(track_dir + "/" + s + "/annotation.html", "w") as file_handle:
        file_handle.write(annotation_html_dic[s])

    with open(track_dir + "/" + s + "/transcriptome.html", "w") as file_handle:
        file_handle.write(transcriptome_html)

    with open(track_dir + "/" + s + "/candidate.html", "w") as file_handle:
        file_handle.write(candidate_html)
    shutil.copyfile("UCSC_track_tools/description_html/" + s + ".jpg",
                    track_dir + "/" + s + "/" + s + ".jpg")
