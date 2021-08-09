#!/usr/bin/python3

'''
This script depends on candidates/make_html.py
'''

import jinja2
import json
import sys
import os
import re
from Bio import SeqIO
from math import log
import difflib
from splinter import Browser
from selenium import webdriver
import subprocess
from glob import glob

from pyteomics import mzml
import matplotlib.pyplot as plt
from spectrum_utils import plot
from spectrum_utils import spectrum


def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


def plot_spectra(mzml_id, peptide, scan_id, mzml_dir, spec_pic_dir, psm_id):
    mzml_file = find(mzml_id, mzml_dir)
    with mzml.read(mzml_file) as reader:
        # auxiliary.print_tree(next(reader))
        for scan in reader:
            if not scan["index"] == int(scan_id) - 1:
                continue
            if "precursorList" not in scan.keys():
                print("no precursor list")
                return
            mz = scan['m/z array']
            intensity = scan['intensity array']
            identifier = scan['index']
            retention_time = float(scan['scanList']['scan'][0]["scan start time"]) * 60.0
            precursor_mz = scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
            precursor_charge = int(scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["charge state"])
            spec = spectrum.MsmsSpectrum(identifier, precursor_mz,
                                         precursor_charge, mz, intensity,
                                         retention_time=retention_time,
                                         peptide=peptide)
            min_mz, max_mz = 100, 1400
            fragment_tol_mass, fragment_tol_mode = 10, 'ppm'
            min_intensity, max_num_peaks = 0.05, 150
            scaling = 'root'
            ion_types = 'aby'
            spec = spec.set_mz_range(min_mz, max_mz)
            spec = spec.remove_precursor_peak(fragment_tol_mass,
                                              fragment_tol_mode)
            spec = spec.filter_intensity(min_intensity, max_num_peaks)
            spec = spec.scale_intensity(scaling)
            spec = spec.annotate_peptide_fragments(fragment_tol_mass,
                                                   fragment_tol_mode,
                                                   ion_types)
            plt.figure()
            plot.spectrum(spec)
            mzml_id = os.path.splitext(os.path.split(mzml_file)[1])[0]
            plt.savefig("{}/{}_{}.svg".format(spec_pic_dir, mzml_id, psm_id))
            plt.close()
            print("print")
        else:
            print("Scan not found")


def system_call(call_list):
    p = subprocess.Popen(call_list, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    p.wait()
    if p.returncode != 0:
        print("error:")
        print(p.communicate()[0].decode('utf-8'))


def render_jinja(file_name, context, template_dir="./"):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir + "/"))
    return env.get_template(file_name).render(context)


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
            if percent > 80:
                del(seqs[j])
            else:
                j += 1
        i += 1
    return len(seqs)


def load_data():
    selection_file = "nov_psm6"

    with open("./parameters.json", "r") as file_handle:
        parameters = json.load(file_handle)
        data_dir = parameters["data_dir"]
        tmp_dir = parameters["tmp_dir"]
        hub_id = parameters["hub_id"]
        session_id = parameters["session_id"]
        chrome_bin = parameters["chrome_bin"]

    url_temp_base = ("https://genome-euro.ucsc.edu/cgi-bin/"
                     "hgTracks?db=hub_{0}_{2}&lastVirtModeType=default&"
                     "lastVirtModeExtraState=&virtModeType=default"
                     "&virtMode=0&nonVirtPosition=&"
                     "position={2}%3A{2}-{2}{1}")
    url_temp = url_temp_base.format(hub_id, session_id, "{}")

    cand_dir = data_dir + "/candidates"
    accu_data_dir = data_dir + "/accumulated_data"
    db_dir = data_dir + "/dbs"
    genome_dir = data_dir + "/genome"
    output_dir = cand_dir + "/html/" + selection_file
    pic_dir = tmp_dir + "/start_anno_pdf_pics"

    if not os.path.exists(pic_dir):
        os.mkdir(pic_dir)

    if not os.path.isdir(output_dir):
        print("Candidate result page does not exist. Build first!")
        sys.exit(1)

    html_template = "start_anno_html/html_templates/start_anno_list.html"

    with open(db_dir + "/orf_to_CDS.json", "r") as file_handle:
        orf_to_CDS = json.load(file_handle)

    frame_dic = SeqIO.index(db_dir + "/SIHUMI_6frame.fasta", "fasta")

    with open(genome_dir + "/uniprot_id_dic.json", "r") as file_handle:
        uniprot_id_dic = json.load(file_handle)

    with open(accu_data_dir + "/earlier_start.json", "r") as file_handle:
        earlier_start = json.load(file_handle)

    with open("../early_starts_new.csv", "r") as file_handle:
        content = [l.split(",") for l in file_handle.read().split("\n")][1:-1]
    category_candidates = {l[0].replace(" ", ""): [l[1].replace(" ", ""), l[2].replace(" ", "")]
                           for l in content}

    chrome_options = set_up_chrome(chrome_bin, pic_dir)

    return (html_template, output_dir, url_temp, orf_to_CDS, frame_dic,
            uniprot_id_dic, earlier_start, pic_dir, chrome_options,
            category_candidates)


def mean_top_psms(psms):
    psm_scores = [psm["e-value"] for psm in psms]
    psm_scores.sort()
    return sum(psm_scores[0:3])/3


def build_context_html(earlier_start, url_temp, orf_to_CDS, frame_dic,
                       uniprot_id_dic, chrome_options, pic_dir,
                       category_candidates):
    context = {}
    context["protein_list"] = [protein for protein, psms in earlier_start.items()
                               if len(psms) > 5]
    as_regex = re.compile("[IL]")
    info_dic = {}
    for protein, psms in earlier_start.items():
        if len(psms) < 6:
            continue
        info = {}
        info["num_psm"] = len(psms)
        info["mean_top_psms"] = log(mean_top_psms(psms), 10) * -1
        info["eval"] = category_candidates[protein][1]
        info["id"] = protein.split("|")[2]
        info["start"] = protein.split("|")[4].split("-")[0]
        info["stop"] = protein.split("|")[4].split("-")[1]
        info["strand"] = protein.split("|")[3]
        info["species"] = protein.split("|")[0]
        info["contig"] = protein.split("|")[1]
        info["ucsc_link"] = url_temp.format("ecoli", "U00096.3",
                                            info["start"], info["stop"])
        info["unique_psms"] = count_unique_psm(psms)
        # get psm which provide evidence for early start
        protein_seq = frame_dic[protein].seq
        start_anno = orf_to_CDS["ecoli"][protein][1] / 3
        early_psms = []
        for psm in psms:
            pep_seq_regex = as_regex.sub("[IL]", psm["pep"])
            starts = [m.start() for m in re.finditer(pep_seq_regex,
                                                     str(protein_seq))]
            if len(starts) != 1:
                print("psm is not ambigous in protein!")
                print(psm)
                print(starts)
            if starts[0] < start_anno:
                early_psms.append(psm)
        info["early_psms"] = early_psms
        # information annotated protein
        annotated_protein = orf_to_CDS["ecoli"][protein][0].split("|")[2]
        uniprot_id = uniprot_id_dic["ecoli"][annotated_protein]
        info["uniprot_id"] = uniprot_id
        info["ncbi_id"] = annotated_protein
        # picture
        # full gene
        get_genome_browser_pic(chrome_options, info["strand"], info["start"],
                               info["stop"], url_temp, info["species"],
                               pic_dir, info["contig"], info["species"] + "_" +
                               info["id"], off_set=0.1)
        # start of gene
        if info["strand"] == "1":
            gene_start = str(int(info["start"]) - 100)
            gene_stop = str(int(info["start"]) + int(start_anno * 3) + 100)
        else:
            gene_start = str(int(info["stop"]) - int(start_anno * 3) - 100)
            gene_stop = str(int(info["stop"]) + 100)
        get_genome_browser_pic(chrome_options, info["strand"], gene_start,
                               gene_stop, url_temp, info["species"], pic_dir,
                               info["contig"], "start" + "_" + info["species"] +
                               "_" + info["id"], off_set=0.1)
        info_dic[protein] = info
    context["protein_info_dic"] = info_dic

    return context


def get_genome_browser_pic(chrome_options, strand, start, stop, url_temp,
                           species, download_dir, contig, name, off_set=0.1):
    print_offset = int((int(stop) - int(start)) * 0.1)
    print_start = int(start) - print_offset
    print_stop = int(stop) + print_offset
    url = url_temp.format(species, contig, print_start, print_stop)
    with Browser("chrome", options=chrome_options) as browser:
        browser.visit(url)
        if strand == "-1":
            element = browser.find_by_xpath('//*[@title="complement bases"]')
            if len(element) != 0:
                element.click()
        pdf_link = browser.find_by_id('pdfLink')
        pdf_url = pdf_link['href']
        browser.visit(pdf_url)
        pdf_link = browser.find_by_text("the current browser graphic in PDF")
        system_call(["wget", "-O", download_dir + "/" + name + ".pdf",
                    pdf_link["href"]])
        if strand == "-1":
            browser.visit(url)
            element = browser.find_by_xpath('//*[@title="complement bases"]')
            if len(element) != 0:
                element.click()


def pdf_to_svg(pic_dir, output_dir):
    output_dir = output_dir + "/pics"
    for pdf_file in glob(pic_dir + "/*pdf"):
        name = pdf_file.split("/")[-1].split(".")[0]
        system_call(["inkscape", "-l", output_dir + "/" + name + ".svg",
                    pic_dir + "/" + name + ".pdf"])


def set_up_chrome(chrome_bin, pic_dir):
    chrome_options = webdriver.ChromeOptions()
    chrome_options.binary_location = chrome_bin
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--no-sandbox")
    prefs = {"download.default_directory": pic_dir,
             "download.directory_upgrade": "true",
             "download.prompt_for_download": "false",
             "disable-popup-blocking": "true"}

    chrome_options.add_experimental_option("prefs", prefs)

    # options that headless works like with window
    chrome_options.add_argument("--disable-infobars")
    chrome_options.add_argument("--window-size=1920,1080")
    chrome_options.add_argument("--start-maximized")
    chrome_options.add_argument("--headless")
    return chrome_options


def print_spectras_PSM(context, mzml_dir, output_dir):
    spec_pic_dir = output_dir + "/pics"
    for protein, info in context["protein_info_dic"].items():
        psms = sorted(info["early_psms"], key=lambda x: x["e-value"])[0:11]
        info["printed_early_psms"] = psms
        for psm in psms:
            plot_spectra(psm["experiment"], psm["pep"], psm["scan"], mzml_dir,
                         spec_pic_dir, str(psm["num"]) + "_" + psm["scan"])


def main():
    (html_template, output_dir, url_temp, orf_to_CDS, frame_dic,
     uniprot_id_dic, earlier_start, pic_dir, chrome_options,
     category_candidates) = load_data()
    context = build_context_html(earlier_start, url_temp, orf_to_CDS,
                                 frame_dic, uniprot_id_dic, chrome_options,
                                 pic_dir, category_candidates)
    pdf_to_svg(pic_dir, output_dir)
    context = print_spectras_PSM(context)

    with open("start_anno_html/context.json", "w") as file_handle:
        json.dump(context, file_handle)

    with open(output_dir + "/start_anno.html", "w") as file_handle:
        file_handle.write(render_jinja(html_template, context))


'''
(html_template, output_dir, url_temp, orf_to_CDS, frame_dic, uniprot_id_dic, earlier_start, pic_dir, chrome_options, category_candidates) = load_data()
context = build_context_html(earlier_start, url_temp, orf_to_CDS, frame_dic, uniprot_id_dic, chrome_options, pic_dir, category_candidates)
pdf_to_svg(pic_dir, output_dir)
context = print_spectras_PSM(context)

with open("start_anno_html/context.json", "w") as file_handle:
    json.dump(context, file_handle)

with open(output_dir + "/start_anno.html", "w") as file_handle:
    file_handle.write(render_jinja(html_template, context))
'''
