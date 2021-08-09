#!/usr/bin/python3

import json
from Bio import SeqIO
from itertools import chain
from copy import deepcopy


def make_orf_dic(frame_dic):
    orf_dic = {}
    orf_list = frame_dic.keys()
    for orf in orf_list:
        species = orf.split('|')[0]
        contig = orf.split("|")[1]
        strand = orf.split("|")[3]
        start = int(orf.split("|")[4].split("-")[0])
        stop = int(orf.split("|")[4].split("-")[1])
        if strand == "1":
            frame = stop % 3
        else:
            frame = start % 3

        if species not in orf_dic:
            orf_dic[species] = {}

        if contig not in orf_dic[species]:
            orf_dic[species][contig] = {}

        if strand not in orf_dic[species][contig]:
            orf_dic[species][contig][strand] = {}

        if frame not in orf_dic[species][contig][strand]:
            orf_dic[species][contig][strand][frame] = []

        orf_dic[species][contig][strand][frame].append(orf)

    sorted_orf_dic = deepcopy(orf_dic)
    for species, species_orf_dic in orf_dic.items():
        for contig, contig_orf_dic in species_orf_dic.items():
            for strand, strand_orf_dic in contig_orf_dic.items():
                for frame, frame_orf_list in strand_orf_dic.items():
                    orf_list = sorted(frame_orf_list, key=lambda protein:
                                      get_start(protein))
                    sorted_orf_dic[species][contig][strand][frame] = orf_list

    return sorted_orf_dic


def get_orf(target_start, target_stop, target_strand, target_contig,
            target_species, orf_dic):
    '''Binary search followed by walking iterativly back to find start of
    intervalls, then walk forward till end of intervalls is found'''
    orf_list = []
    if target_strand == "-1":
        target_frame = target_start % 3
    else:
        target_frame = target_stop % 3
    all_orfs = orf_dic[target_species][target_contig][target_strand][target_frame]
    smallest_i = 0
    highest_i = len(all_orfs)
    i = int(len(all_orfs) / 2)
    binary_finished = False
    start_found = False
    found = False
    while True:
        # with joined features it is possible that the iterator over or under
        # reaches the list
        if i < 0 or i > len(all_orfs):
            return
        orf = all_orfs[i]
        start = int(orf.split("|")[4].split("-")[0])
        stop = int(orf.split("|")[4].split("-")[1])

        # first binary
        # first condition checks if binary seach is finished, but start and
        # first orfs are collected
        if binary_finished and not found and not start_found:
            pass
        elif target_stop < start:
            # if found all orfs are collected finished!
            if found:
                break
            highest_i = i
            i = int((smallest_i + i) / 2)
            continue
        # or intervall is lower
        elif stop < target_start:
            if found:
                break
            smallest_i = i
            i = int((highest_i + i) / 2)
            continue

        # after jumps in intervall next find start of intervall
        if not binary_finished:
            binary_finished = True
            if i == 0:
                start_found = True
            else:
                i -= 1

        if not start_found:
            if stop < target_start:
                start_found = True
                i += 1
            else:
                i -= 1
            continue

        # after start of intervall found collect all orfs

        orf_list.append(orf)
        found = True
        i += 1
    return orf_list


def get_sub_protein(CDS, orf):
    CDS_start = int(CDS.split("|")[4].split("-")[0])
    CDS_stop = int(CDS.split("|")[4].split("-")[1])
    orf_start = int(orf.split("|")[4].split("-")[0])
    orf_stop = int(orf.split("|")[4].split("-")[1])
    orf_strand = orf.split("|")[3]
    CDS_strand = CDS.split("|")[3]

    if orf_strand != CDS_strand:
        print("Error CDS: {} and orf: {} not on same strand".format(CDS, orf))
        raise ValueError
    if not ((orf_start <= CDS_stop) and (CDS_start <= orf_stop)):
        print("Error CDS: {} and orf: {} not overlapping".format(CDS, orf))
        raise ValueError
    if orf_strand == "1":
        sub_start = CDS_start - orf_start
    else:
        sub_start = orf_stop - CDS_stop
    return sub_start


def get_start(protein):
    return int(protein.split("|")[4].split("-")[0])


def check_CDS(proteom_dic, frame_dic, orf_dic):
    CDS_to_orf = {}
    orf_to_CDS = {}
    CDS_intervall_dic = {}
    anno_orf = []
    print("checking each CDS")

    for CDS, CDS_seq in proteom_dic.items():
        CDS_species = CDS.split("|")[0]
        CDS_contig = CDS.split("|")[1]
        CDS_start = int(CDS.split("|")[4].split("-")[0])
        CDS_stop = int(CDS.split("|")[4].split("-")[1])
        CDS_strand = CDS.split("|")[3]
        if CDS_species not in CDS_to_orf:
            CDS_to_orf[CDS_species] = {}

        if CDS_species not in orf_to_CDS:
            orf_to_CDS[CDS_species] = {}

        if CDS_species not in CDS_intervall_dic:
            CDS_intervall_dic[CDS_species] = {}

        if CDS_contig not in CDS_intervall_dic[CDS_species]:
            CDS_intervall_dic[CDS_species][CDS_contig] = {}

        if CDS_strand not in CDS_intervall_dic[CDS_species][CDS_contig]:
            CDS_intervall_dic[CDS_species][CDS_contig][CDS_strand] = []

        orf_proteins = get_orf(CDS_start, CDS_stop, CDS_strand, CDS_contig,
                               CDS_species, orf_dic)

        if orf_proteins is None:
            CDS_to_orf[CDS_species][CDS] = [orf_proteins, "annotation error"]
            continue
        elif (len(orf_proteins) != 1):
            CDS_to_orf[CDS_species][CDS] = [orf_proteins, "annotation error"]
            CDS_intervall_dic[CDS_species][CDS_contig][CDS_strand].append([CDS_start,
                                                                           CDS_stop])
            anno_orf.append(orf_proteins)
            for orf in orf_proteins:
                orf_to_CDS[CDS_species][orf] = [CDS, "annotation error",
                                                "intragenic"]
            continue

        start_anno = get_sub_protein(CDS, orf_proteins[0])
        start_anno_prot = int(start_anno / 3)

        # check exception
        CDS_sub_seq = str(CDS_seq[1:].seq)
        orf_sub_seq = str(frame_dic[orf_proteins[0]].seq[1 + start_anno_prot:])
        if orf_sub_seq != CDS_sub_seq:
            CDS_to_orf[CDS_species][CDS] = [orf_proteins, "annotation error"]
            orf_to_CDS[CDS_species][orf_proteins[0]] = [CDS, "annotation error",
                                                        "intragenic"]
            CDS_intervall_dic[CDS_species][CDS_contig][CDS_strand].append([CDS_start,
                                                                           CDS_stop])
            anno_orf.append(orf_proteins)
            continue

        # normal orf to CDS relation goes here
        CDS_to_orf[CDS_species][CDS] = [orf_proteins, start_anno]
        orf_to_CDS[CDS_species][orf_proteins[0]] = [CDS, start_anno, "intragenic"]
        CDS_intervall_dic[CDS_species][CDS_contig][CDS_strand].append([CDS_start,
                                                                       CDS_stop])
        anno_orf.append(orf_proteins)

    anno_orf = list(chain.from_iterable(anno_orf))

    sorted_CDS_intervall_dic = deepcopy(CDS_intervall_dic)
    for species, species_CDS_intervall_dic in CDS_intervall_dic.items():
        for contig, contig_CDS_intervall_dic in species_CDS_intervall_dic.items():
            for strand, intervall_list in contig_CDS_intervall_dic.items():
                orf_list = sorted(intervall_list)
                sorted_CDS_intervall_dic[species][contig][strand] = orf_list

    return orf_to_CDS, CDS_to_orf, anno_orf, sorted_CDS_intervall_dic


def check_orf(frame_dic, orf_to_CDS, CDS_to_orf, anno_orf, CDS_intervall_dic):
    print("checking orf")
    for orf in frame_dic.keys():
        if orf in anno_orf:
            continue
        species = orf.split("|")[0]
        contig = orf.split("|")[1]
        strand = orf.split("|")[3]
        start = int(orf.split("|")[4].split("-")[0])
        stop = int(orf.split("|")[4].split("-")[1])
        genomic_postion = "intergenic"

        if species not in orf_to_CDS:
            orf_to_CDS[species] = {}

        # x1 <= y2 && y1 <= x2
        for interval in CDS_intervall_dic[species][contig][strand]:
            if interval[1] < start:
                continue
            if interval[0] > stop:
                break
            genomic_postion = "same strand intragenic"
            break
        if genomic_postion == "same strand intragenic":
            orf_to_CDS[species][orf] = [None, None, genomic_postion]
            continue

        for interval in CDS_intervall_dic[species][contig][str(int(strand) * -1)]:
            if interval[1] < start:
                continue
            if interval[0] > stop:
                break
            genomic_postion = "opposite strand intragenic"
            break

        orf_to_CDS[species][orf] = [None, None, genomic_postion]
    return orf_to_CDS


def main():
    print("loading data")
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    db_dir = data_dir + "/dbs"
    frame_path = db_dir + "/SIHUMI_6frame.fasta"
    frame_dic = SeqIO.index(frame_path, "fasta")

    proteom_path = db_dir + "/SIHUMI_proteom.fasta"
    proteom_dic = SeqIO.index(proteom_path, "fasta")

    orf_dic = make_orf_dic(frame_dic)
    orf_to_CDS, CDS_to_orf, anno_orf, CDS_intervall_dic = check_CDS(proteom_dic,
                                                                    frame_dic,
                                                                    orf_dic)

    orf_to_CDS = check_orf(frame_dic, orf_to_CDS, CDS_to_orf, anno_orf,
                           CDS_intervall_dic)

    with open(db_dir + "/CDS_to_orf.json", "w") as file_handle:
        json.dump(CDS_to_orf, file_handle)

    with open(db_dir + "/orf_to_CDS.json", "w") as file_handle:
        json.dump(orf_to_CDS, file_handle)


if __name__ == "__main__":
    main()


"""
"""
