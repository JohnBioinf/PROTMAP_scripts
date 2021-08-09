#!/usr/bin/python3

import json
import os
from html.parser import HTMLParser
from Bio import SeqIO
from urllib.request import urlopen
from time import sleep
import re
import requests


def get_curent_uniprotID(uniprot_id):
    if type(uniprot_id) == list:
        results = []
        for u in uniprot_id:
            results.append(get_curent_uniprotID(u))

        return results

    if uniprot_id == "-":
        return "-"

    tries = 1
    err = ""
    pars = HTMLParser()
    while tries < 6:
        try:
            response = urlopen("https://www.uniprot.org/uniprot/" + uniprot_id)
            uniprot_html = pars.unescape(response.read().decode("utf8"))
            err = None
        except Exception as error:
            if str(err) == "HTTP Error 404: ":
                print("not a uniport id: " + uniprot_id)
                return "-"
            sleep(2)
            err = error
            tries = tries + 1
        if err is None:
            break
    else:
        print("https://www.uniprot.org/uniprot/" + uniprot_id)
        print(err)
        raise Exception(err)

    if "This entry is obsolete" in uniprot_html:
        return None
    else:
        return response.url.split("/")[-1]


def get_validation_uniprotID(uniprot_id):
    if type(uniprot_id) == list:
        results = []
        for u in uniprot_id:
            results.append(get_validation_uniprotID(u))
        return results
    validation_dic = {"5": "Experimental evidence at protein level",
                      "4": "Experimental evidence at transcript level",
                      "3": "Protein inferred from homology",
                      "2": "Protein predicted",
                      "1": "Protein uncertain"}
    tries = 1
    err = ""
    while tries < 6:
        try:
            response = urlopen("https://www.uniprot.org/uniprot/" + uniprot_id)
            uniprot_html = response.read().decode("utf8")
            err = None
        except Exception as error:
            if str(error) == "HTTP Error 404: ":
                print("not a uniport id: " + uniprot_id)
                return "-"
            sleep(2)
            err = error
            tries = tries + 1
        if err is None:
            break
    else:
        print("https://www.uniprot.org/uniprot/" + uniprot_id)
        raise Exception(err)
    if "This entry is obsolete" in uniprot_html:
        return "-"
    else:
        line = [l for l in uniprot_html.split("\n") if "<p>Annotation score:" in l][0]
        validation_value = re.split("score:| ", re.search("<p>Annotation score:.*</p>", line).group(0))[2]
        return validation_value + " " + validation_dic[validation_value]


def db2db(input_ids, input_type="genesymbol", output_type="geneid", taxid=None):
    output_id_list = ["Affy ID", "Agilent ID", "Allergome Code",
                      "ApiDB_CryptoDB ID", "Biocarta Pathway Name", "BioCyc ID",
                      "CCDS ID", "Chromosomal Location", "CleanEx ID",
                      "CodeLink ID", "COSMIC ID", "CPDB Protein Interactor",
                      "CTD Disease Info", "CTD Disease Name", "CYGD ID",
                      "dbSNP ID", "dictyBase ID", "DIP ID", "DisProt ID",
                      "DrugBank Drug ID", "DrugBank Drug Info",
                      "DrugBank Drug Name", "EC Number", "EchoBASE ID",
                      "EcoGene ID", "Ensembl Biotype", "Ensembl Gene ID",
                      "Ensembl Gene Info", "Ensembl Protein ID",
                      "Ensembl Transcript ID", "FlyBase Gene ID",
                      "FlyBase Protein ID", "FlyBase Transcript ID",
                      "GAD Disease Info", "GAD Disease Name",
                      "GenBank Nucleotide Accession", "GenBank Nucleotide GI",
                      "GenBank Protein Accession", "GenBank Protein GI",
                      "Gene Info", "Gene Symbol", "Gene Symbol and Synonyms",
                      "Gene Symbol ORF", "Gene Synonyms", "GeneFarm ID",
                      "GO - Biological Process", "GO - Cellular Component",
                      "GO - Molecular Function", "GO ID", "GSEA Standard Name",
                      "H-Inv Locus ID", "HAMAP ID", "HGNC ID",
                      "HMDB Metabolite", "Homolog - All Ens Gene ID",
                      "Homolog - All Ens Protein ID", "Homolog - All Gene ID",
                      "Homolog - Human Ens Gene ID",
                      "Homolog - Human Ens Protein ID",
                      "Homolog - Human Gene ID", "Homolog - Mouse Ens Gene ID",
                      "Homolog - Mouse Ens Protein ID",
                      "Homolog - Mouse Gene ID", "Homolog - Rat Ens Gene ID",
                      "Homolog - Rat Ens Protein ID", "Homolog - Rat Gene ID",
                      "HomoloGene ID", "HPA ID", "HPRD ID",
                      "HPRD Protein Complex", "HPRD Protein Interactor",
                      "Illumina ID", "IMGT/GENE-DB ID", "InterPro ID", "IPI ID",
                      "KEGG Disease ID", "KEGG Gene ID", "KEGG Orthology ID",
                      "KEGG Pathway ID", "KEGG Pathway Info",
                      "KEGG Pathway Title", "LegioList ID", "Leproma ID",
                      "Locus Tag", "MaizeGDB ID", "MEROPS ID",
                      "MGC(ZGC/XGC) ID", "MGC(ZGC/XGC) Image ID",
                      "MGC(ZGC/XGC) Info", "MGI ID", "MIM ID", "MIM Info",
                      "miRBase ID", "NCIPID Pathway Name",
                      "NCIPID Protein Complex", "NCIPID Protein Interactor",
                      "NCIPID PTM", "Orphanet ID", "PANTHER ID",
                      "Paralog - Ens Gene ID", "PBR ID", "PDB ID",
                      "PeroxiBase ID", "Pfam ID", "PharmGKB Drug Info",
                      "PharmGKB Gene ID", "PIR ID", "PIRSF ID", "PptaseDB ID",
                      "PRINTS ID", "ProDom ID", "PROSITE ID", "PseudoCAP ID",
                      "PubMed ID", "Reactome ID", "Reactome Pathway Name",
                      "REBASE ID", "RefSeq Genomic Accession",
                      "RefSeq Genomic GI", "RefSeq mRNA Accession",
                      "RefSeq ncRNA Accession", "RefSeq Nucleotide GI",
                      "RefSeq Protein Accession", "RefSeq Protein GI",
                      "Rfam ID", "RGD ID", "SGD ID", "SMART ID",
                      "STRING Protein Interactor", "TAIR ID", "Taxon ID",
                      "TCDB ID", "TIGRFAMs ID", "TubercuList ID", "UCSC ID",
                      "UniGene ID", "UniProt Accession", "UniProt Entry Name",
                      "UniProt Info", "UniProt Protein Name", "UniSTS ID",
                      "VectorBase Gene ID", "VEGA Gene ID", "VEGA Protein ID",
                      "VEGA Transcript ID", "WormBase Gene ID",
                      "WormPep Protein ID", "XenBase Gene ID", "ZFIN ID"]

    input_id_list = ["Affy GeneChip Array", "Affy ID",
                     "Affy Transcript Cluster ID", "Agilent ID",
                     "Biocarta Pathway Name", "CodeLink ID", "dbSNP ID",
                     "DrugBank Drug ID", "DrugBank Drug Name", "EC Number",
                     "Ensembl Gene ID", "Ensembl Protein ID",
                     "Ensembl Transcript ID", "EST Accession",
                     "FlyBase Gene ID", "GenBank Nucleotide Accession",
                     "GenBank Protein Accession", "Gene ID", "Gene Symbol",
                     "Gene Symbol and Synonyms", "Gene Symbol ORF", "GI Number",
                     "GO ID", "GSEA Standard Name", "H-Inv Locus ID",
                     "H-Inv Protein ID", "H-Inv Transcript ID", "HGNC ID",
                     "HMDB Metabolite", "HomoloGene ID", "Illumina ID",
                     "InterPro ID", "IPI ID", "KEGG Compound ID",
                     "KEGG Compound Name", "KEGG Disease ID", "KEGG Drug ID",
                     "KEGG Drug Name", "KEGG Gene ID", "KEGG Pathway ID",
                     "MaizeGDB ID", "MGI ID", "MIM ID", "miRBase ID",
                     "miRBase Mature miRNA Acc", "NCIPID Pathway Name",
                     "Organism Scientific Name", "PDB ID", "Pfam ID",
                     "PharmGKB ID", "PubChem ID", "Reactome Pathway Name",
                     "RefSeq Genomic Accession", "RefSeq mRNA Accession",
                     "RefSeq Protein Accession", "SGD ID", "TAIR ID",
                     "Taxon ID", "UniGene ID", "UniProt Accession",
                     "UniProt Entry Name", "UniProt Protein Name", "UniSTS ID"]
    output_type_list = [i.lower().replace(" ", "") for i in output_id_list]
    input_type_list = [i.lower().replace(" ", "") for i in input_id_list]
    input_type = input_type.lower().replace(" ", "")
    output_type = output_type.lower().replace(" ", "")
    if output_type not in output_type_list:
        raise Exception("Wrong output type: {}.".format(output_type))
    result_key = output_id_list[output_type_list.index(output_type)]
    if input_type not in input_type_list:
        raise Exception("Wrong input type {}.".format(input_type))
    max_input = 500
    # clean version nummer Possible problem
    # version number start
    if type(input_ids) is list:
        has_version_number = any([i.split(".")[-1].isdigit() and
                                  len(i.split(".")) > 1
                                  for i in input_ids])
    else:
        has_version_number = input_ids.isdigit()

    if input_type == "refseqproteinaccession" and has_version_number:
        if type(input_ids) is list:
            input_ids = [".".join(i.split(".")[0:-1])
                         for i in input_ids
                         if i.split(".")[-1].isdigit() and
                         len(i.split(".")) > 1]
        else:
            input_ids = ".".join(input_ids.split(".")[0:-1])
    elif has_version_number:
        if type(input_ids) is list:
            version_ids = [i for i in input_ids
                           if i.split(".")[-1].isdigit() and
                           len(i.split(".")) > 1]
            print("Ids have version number. Possible problem.")
            print(input_type)
            print("\n".join(version_ids))
        else:
            print("Id " + input_ids + " has version number. Possible problem.")

    # version number end

    if type(input_ids) is list:
        if len(input_ids) > max_input:
            results = []
            for i in range(0, len(input_ids), max_input):
                results = results + db2db(input_ids[i:min(i + max_input, len(input_ids))], input_type, output_type, taxid)
            return results
        else:
            input_ids = ",".join(input_ids)

    if taxid is not None:
        taxid = "&taxonId=" + taxid
    else:
        taxid = ""

    base_url = "https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input={}&inputValues={}&outputs={}"
    url = base_url.format(input_type, input_ids, output_type) + taxid

    result = requests.get(url)
    result = result.json()
    return [r[result_key].split("//") for r in result]


def get_uniprot_new_seq(genome, uniprot_dic):
    print("\nGet uniprot for new sequenced species")
    species_list = ["lacto", "anaero", "blautia", "bact", "ery", "clostri"]
    refseq_dic = {}
    for species in species_list:
        print(species)
        refseq_dic[species] = {}
        for contig_name, contig in genome[species].items():
            for feature in contig.features:
                if feature.type == "CDS":
                    inferences = feature.qualifiers["inference"]
                    if "protein_id" not in feature.qualifiers.keys():
                        if "locus_tag" in feature.qualifiers.keys():
                            protein_id = feature.qualifiers["locus_tag"][0]
                    else:
                        protein_id = feature.qualifiers["protein_id"][0]
                    if any(["RefSeq" in i for i in inferences]):
                        inferences = [i for i in inferences if "RefSeq" in i]
                        if len(inferences) > 1:
                            print("Unknown inference combination")
                            print(species)
                            print(feature)
                            raise
                        RefSeq_id = inferences[0].split("COORDINATES: similar to AA sequence:")[-1].replace("RefSeq:", "").split(",")[0]
                        refseq_dic[species][protein_id] = RefSeq_id
                    else:
                        if len(inferences) > 1:
                            print("Unknown inference combination")
                            print(species)
                            print(feature)
                            raise
                        if "HMM" in inferences[0]:
                            refseq_dic[species][protein_id] = "-"
                        elif "ab initio prediction" in inferences[0]:
                            refseq_dic[species][protein_id] = "-"
                        else:
                            print("Unknown inference")
                            print(species)
                            print(feature)
                            raise

    print("\nCheck RefSeq for Uniprot")
    uniprot_id_dic = {}
    for species, species_refseq_dic in refseq_dic.items():
        print(species)
        uniprot_id_dic[species] = {}
        RefSeq_ids = [i for prot_id, i in species_refseq_dic.items() if i != "-"]
        uniprot_id_list = db2db(RefSeq_ids,
                                input_type="RefSeq Protein Accession",
                                output_type="UniProt Accession")
        uniprot_refseq_dic = {}
        uniprot_refseq_dic["-"] = "-"
        for i in range(0, len(uniprot_id_list)):
            uniprot_refseq_dic[RefSeq_ids[i]] = uniprot_id_list[i]

        for protein_id, RefSeq_id in refseq_dic[species].items():
            uniprot_id_dic[species][protein_id] = uniprot_refseq_dic[RefSeq_id]

    return uniprot_id_dic


def get_uniprot_ecoli(genome, uniprot_id_dic):
    species = "ecoli"
    uniprot_id_dic[species] = {}
    kegg_dic = {}
    for contig_name, contig in genome[species].items():
        for feature in contig.features:
            if feature.type == "CDS":
                if "protein_id" not in feature.qualifiers.keys():
                    if "locus_tag" in feature.qualifiers.keys():
                        protein_id = feature.qualifiers["locus_tag"][0]
                else:
                    protein_id = feature.qualifiers["protein_id"][0]
                for ref in feature.qualifiers["db_xref"]:
                    if "UniProtKB/Swiss-Prot" in ref:
                        uniprot_id_dic[species][protein_id] = ref.split(":")[1]
                        break
                else:
                    kegg_dic[protein_id] = "eco:" + feature.qualifiers["locus_tag"][0]

    kegg_ids = [i for prot_id, i in kegg_dic.items()]
    uniprot_id_list = db2db(kegg_ids,
                            input_type="KEGG Gene ID",
                            output_type="UniProt Accession")
    uniprot_kegg_dic = {}
    for i in range(0, len(uniprot_id_list)):
        uniprot_kegg_dic[kegg_ids[i]] = uniprot_id_list[i]

    for protein_id, kegg_id in kegg_dic.items():
        uniprot_id_dic[species][protein_id] = uniprot_kegg_dic[kegg_id]

    return uniprot_id_dic


def get_uniprot_bifi(genome, uniprot_id_dic):
    species = "bifi"
    uniprot_id_dic[species] = {}
    genbank_dic = {}
    for contig_name, contig in genome[species].items():
        for feature in contig.features:
            if feature.type == "CDS":
                if "protein_id" not in feature.qualifiers.keys():
                    if "locus_tag" in feature.qualifiers.keys():
                        protein_id = feature.qualifiers["locus_tag"][0]
                else:
                    protein_id = feature.qualifiers["protein_id"][0]
                genbank_dic[protein_id] = feature.qualifiers["protein_id"][0].split(".")[0]

    genbank_ids = [i for prot_id, i in genbank_dic.items()]
    uniprot_id_list = db2db(genbank_ids,
                            input_type="GenBank Protein Accession",
                            output_type="UniProt Accession")
    uniprot_genbank_dic = {}
    for i in range(0, len(uniprot_id_list)):
        uniprot_genbank_dic[genbank_ids[i]] = uniprot_id_list[i]

    for protein_id, genbank_id in genbank_dic.items():
        uniprot_id_dic[species][protein_id] = uniprot_genbank_dic[genbank_id]

    return uniprot_id_dic


def clean_uniprot(uniprot_id_dic):
    uniprot_id_dic_clean = {}
    for species, species_uniprot_id_dic in uniprot_id_dic.items():
        uniprot_id_dic_clean[species] = {}
        for protein_id, uniprot_id_list in species_uniprot_id_dic.items():
            if uniprot_id_list == "-":
                uniprot_id_dic_clean[species][protein_id] = ["-"]
            uniprot_id_list_clean = get_curent_uniprotID(uniprot_id_list)

            if type(uniprot_id_list) is str:
                uniprot_id_list_clean = [uniprot_id_list_clean]

            uniprot_id_list_clean = [i for i in uniprot_id_list_clean if i is not None]

            uniprot_id_dic_clean[species][protein_id] = list(set(uniprot_id_list_clean))

    return uniprot_id_dic_clean


def get_validation(uniprot_id_dic):
    validation_dic = {}
    for species, species_uniprot_id_dic in uniprot_id_dic.items():
        validation_dic[species] = {}
        for protein_id, uniprot_id_list in species_uniprot_id_dic.items():
            if type(uniprot_id_list) is list and len(uniprot_id_list) > 1:
                if type(uniprot_id_list[0]) == str and len(uniprot_id_list[0]) == 1:
                    print(uniprot_id_list)
                    print(protein_id)
                    print(species)
                    exit()
            validation_levels = []
            for uniprot_id in uniprot_id_list:
                if (uniprot_id == "-"):
                    validation_dic[species][protein_id] = "0 Protein hypothetical"
                    continue
                validation_levels.append(get_validation_uniprotID(uniprot_id))
            if len(validation_levels) == 0:
                validation_dic[species][protein_id] = "0 Protein hypothetical"
                continue
            max_validation_level = sorted(validation_levels, key=lambda v: str(v.split(" ")[0]))[-1]
            validation_dic[species][protein_id] = max_validation_level

    return validation_dic


def main():
    print("loading data")
    with open("./parameters.json", "r") as file_handle:
        param = json.load(file_handle)
        data_dir = param['data_dir']

    genome_dir = data_dir + "/genome"
    genome = {}
    for genbank_file in os.listdir(genome_dir):
        if genbank_file.endswith(".gbk"):
            species = genbank_file.split(".")[0]
            genome_path = genome_dir + "/" + genbank_file
            genome[species] = SeqIO.to_dict(SeqIO.parse(genome_path, "genbank"))
    uniprot_id_dic = {}
    uniprot_id_dic = get_uniprot_new_seq(genome, uniprot_id_dic)
    uniprot_id_dic = get_uniprot_ecoli(genome, uniprot_id_dic)
    uniprot_id_dic = get_uniprot_bifi(genome, uniprot_id_dic)

    with open(genome_dir + "/uniprot_id_dic_dirty.json", "w") as file_handle:
        json.dump(uniprot_id_dic, file_handle)

    with open(genome_dir + "/uniprot_id_dic_dirty.json", "r") as file_handle:
        uniprot_id_dic = json.load(file_handle)

    print("clean uniprot")

    uniprot_id_dic = clean_uniprot(uniprot_id_dic)
    with open(genome_dir + "/uniprot_id_dic.json", "w") as file_handle:
        json.dump(uniprot_id_dic, file_handle)

    with open(genome_dir + "/uniprot_id_dic.json", "r") as file_handle:
        uniprot_id_dic = json.load(file_handle)

    print("get validation")
    validation_dic = get_validation(uniprot_id_dic)

    with open(genome_dir + "/validation_dic.json", "w") as file_handle:
        json.dump(validation_dic, file_handle)

    print("finished")
    return


if __name__ == "__main__":
    main()
