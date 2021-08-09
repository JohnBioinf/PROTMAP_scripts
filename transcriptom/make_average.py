#!/usr/bin/python3

import json
from glob import glob


def load_wigs(wiggle_files):
    sum_wiggle = {}
    with open(wiggle_files[0], "r") as file_handler:
        for line in file_handler:
            line = line[:-1].split(" ")
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            value = float(line[3])
            if stop - start == 1:
                key = chrom + "__" + str(start)
                sum_wiggle[key] = value
            else:
                for i in range(0, stop - start - 1):
                    key = chrom + "__" + str(start + i)
                    sum_wiggle[key] = value

    for wiggle_file in wiggle_files[1:]:
        with open(wiggle_file, "r") as file_handler:
            for line in file_handler:
                line = line[:-1].split(" ")
                chrom = line[0]
                start = int(line[1])
                stop = int(line[2])
                value = float(line[3])
                if stop - start == 1:
                    key = chrom + "__" + str(start)
                    '''
                    if key == "U00096.3__185":
                        print(value)
                        print(sum_wiggle[key])
                        print(line)
                    '''
                    if key in sum_wiggle:
                        sum_wiggle[key] += value
                    else:
                        sum_wiggle[key] = value
                else:
                    for i in range(0, stop - start - 1):
                        key = chrom + "__" + str(start + i)
                        if key in sum_wiggle:
                            sum_wiggle[key] += value
                        else:
                            sum_wiggle[key] = value

    return sum_wiggle


def main():
    with open("./parameters.json", "r") as file_handle:
        data_dir = json.load(file_handle)['data_dir']
    wiggle_files = glob(data_dir + "/transcriptom/*.fastq")
    wiggle_files = [f.replace("fastq", "wig") for f in wiggle_files]

    average_file = data_dir + "/transcriptom/average.wig"

    sum_wiggle = load_wigs(wiggle_files)

    print("finished loading wigs")

    result = []
    for key, value in sum_wiggle.items():
        chrom = key.split("__")[0]
        start = int(key.split("__")[1])
        result.append([chrom, str(start), str(start + 1),
                       str(int(value / len(wiggle_files)))])

    print("finished aggregating results")

    result.sort(key=lambda x: (x[0], int(x[1])))

    with open(average_file, "w") as file_handler:
        file_handler.write("\n".join([" ".join(r) for r in result]) + "\n")

    print("finished writing results")


if __name__ == "__main__":
    main()


"""
"""
