#!/usr/bin/python3

from splinter import Browser
from selenium import webdriver
from os import system


def download_pic(hub_id, species, contig, start, stop, session_id, strand,
                 download_dir, name):
    url_temp = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hub_{}_{}&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position={}%3A{}-{}{}"
    url = url_temp.format(hub_id, species, contig, start, stop, session_id)

    # set download_dir
    chrome_options = webdriver.ChromeOptions()
    prefs = {"download.default_directory": download_dir,
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
        browser.visit(url)
        if strand == "-1":
            browser.find_by_xpath('//*[@title="complement bases"]').click()
        pdf_link = browser.find_by_id('pdfLink')
        pdf_url = pdf_link['href']
        browser.visit(pdf_url)
        pdf_link = browser.find_by_text("the current browser graphic in PDF")
        system("wget -O {}/{}.pdf {}".format(download_dir, name,
                                             pdf_link["href"]))
        if strand == "-1":
            browser.visit(url)
            browser.find_by_xpath('//*[@title="complement bases"]').click()
