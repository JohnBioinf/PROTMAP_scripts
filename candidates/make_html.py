#!/usr/bin/python3

import jinja2
import json
from shutil import copyfile
import sys


def render_jinja(file_name, context, template_dir="./"):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir + "/"))
    return env.get_template(file_name).render(context)


def main():
    selection_cut_off = int(sys.argv[1])
    selection_file = "nov_psm" + str(selection_cut_off)

    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']

    main_page_template = "candidates/html_templates/main_page_template.html"
    overview_template = "candidates/html_templates/overview.html"
    java_script_path = "candidates/html_templates/fun.js"
    css_path = "candidates/html_templates/style.css"

    cand_dir = data_dir + "/candidates"
    output_dir = cand_dir + "/html/" + selection_file

    copyfile(java_script_path, output_dir + "/fun.js")
    copyfile(css_path, output_dir + "/style.css")

    with open(cand_dir + "/" + selection_file + "_info_dic.json", "r") as file_handler:
        info_dic = json.load(file_handler)

    context = {}
    context["num_cand"] = len(info_dic)
    context["num_true_cand"] = len([1 for k, v in info_dic.items()
                                    if "True" in v["classification"]])
    context["num_true_novel_cand"] = len([1 for k, v in info_dic.items()
                                          if "True" in v["classification"] and
                                          v["blast_category"] in ["Hypothetical", "Novel"]])
    context["num_true_novel_no_blast_cand"] = len([1 for k, v in info_dic.items()
                                                   if "True" in v["classification"] and
                                                   v["blast_category"] == "Novel"])
    context["num_novel"] = len([1 for k, v in info_dic.items()
                                if v["blast_category"] in ["Hypothetical", "Novel"]])
    context["num_hypo"] = len([1 for k, v in info_dic.items()
                               if v["blast_category"] == "Hypothetical"])
    context["num_high_psm"] = len([1 for k, v in info_dic.items() if
                                   v["mean_top_psms"] > 3.5])

    with open(output_dir + "/overview.html", "w") as file_handler:
        file_handler.write(render_jinja(overview_template, context))

    context = {}
    context["candidate_dic"] = info_dic
    candidate_list = sorted([[int(v["name"].split("_")[-1]), k] for k, v in info_dic.items()])
    context["candidate_list"] = [c[1] for c in candidate_list]

    with open(output_dir + "/cand_list.html", "w") as file_handler:
        file_handler.write(render_jinja(main_page_template, context))


if __name__ == "__main__":
    main()


'''
selection_cut_off = 6
selection_file = "nov_psm" + str(selection_cut_off)

with open("./parameters.json", "r") as file_handle:
    data_dir = json.load(file_handle)['data_dir']

main_page_template = "candidates/html_templates/main_page_template.html"
java_script_path = "candidates/html_templates/fun.js"
css_path = "candidates/html_templates/style.css"

cand_dir = data_dir + "/candidates"
output_dir = cand_dir + "/html/" + selection_file

with open(cand_dir + "/" + selection_file + "_info_dic.json", "r") as file_handler:
    info_dic = json.load(file_handler)

context = {}
context["candidate_dic"] = info_dic
candidate_list = sorted([[int(v["name"].split("_")[-1]), k] for k, v in info_dic.items()])
context["candidate_list"] = [c[1] for c in candidate_list]
'''
