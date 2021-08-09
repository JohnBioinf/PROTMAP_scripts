#!/usr/bin/python3

import json
import sys
from splinter import Browser
from selenium import webdriver
import os


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
        os.system("wget -O {}/{}.pdf {}".format(download_dir, name,
                                                pdf_link["href"]))
        if strand == "-1":
            browser.visit(url)
            element = browser.find_by_xpath('//*[@title="complement bases"]')
            if len(element) != 0:
                element.click()


def main():
    selection_cut_off = int(sys.argv[1])
    selection_file = "nov_psm" + str(selection_cut_off)

    with open("./parameters.json", "r") as file_handle:
        parameters = json.load(file_handle)
        data_dir = parameters["data_dir"]
        tmp_dir = parameters["tmp_dir"]
        hub_id = parameters["hub_id"]
        session_id = parameters["session_id"]
        chrome_bin = parameters["chrome_bin"]

    cand_dir = data_dir + "/candidates"
    pic_dir = tmp_dir + "/candidates_pdf_pics"

    if not os.path.exists(pic_dir):
        os.mkdir(pic_dir)

    with open(cand_dir + "/" + selection_file + "_info_dic.json", "r") as file_handler:
        info_dic = json.load(file_handler)

    with open(cand_dir + "/" + selection_file + "_list.json", "r") as file_handler:
        cand_list = json.load(file_handler)

    url_temp_base = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hub_{0}_{2}&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position={2}%3A{2}-{2}{1}"
    url_temp = url_temp_base.format(hub_id, session_id, "{}")
    # set download_dir
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

    for candidate in cand_list:
        info = info_dic[candidate]
        name = info["name"] + "_gene"
        print(name)
        get_genome_browser_pic(chrome_options, info["strand"], info["start"],
                               info["stop"], url_temp, info["species"], pic_dir,
                               info["contig"], name, off_set=0.1)

        name = info["name"] + "_context"
        print(name)
        get_genome_browser_pic(chrome_options, info["strand"],
                               info["start_genomic_context"],
                               info["stop_genomic_context"], url_temp,
                               info["species"], pic_dir, info["contig"], name,
                               off_set=0.1)


if __name__ == "__main__":
    main()


'''
import json
import sys
from splinter import Browser
from selenium import webdriver
import os

chrome_options = webdriver.ChromeOptions()
chrome_options.binary_location = "/root/bin/chrome-linux/chrome"
chrome_options.add_argument("--disable-dev-shm-usage")
chrome_options.add_argument("--no-sandbox")
prefs = {"download.default_directory": ".",
         "download.directory_upgrade": "true",
         "download.prompt_for_download": "false",
         "disable-popup-blocking": "true"}

chrome_options.add_experimental_option("prefs", prefs)

# options that headless works like with window
chrome_options.add_argument("--disable-infobars")
chrome_options.add_argument("--window-size=1920,1080")
chrome_options.add_argument("--start-maximized")
chrome_options.add_argument("--headless")

with Browser("chrome", options=chrome_options) as browser:
    browser.visit("https://www.google.de")

'''
