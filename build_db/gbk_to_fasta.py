#!/usr/bin/python

from Bio import SeqIO
import sys
import os.path

if len(sys.argv) == 1:
    print("No Arguments")
    sys.exit(1)

input_file = sys.argv[1]
if len(sys.argv) == 3:
    output_file = sys.argv[2]
else:
    output_file = ("/".join(input_file.split("/")[0:-1]) + "/" +
                   input_file.split("/")[-1].split(".")[0] + ".fasta")

if not os.path.isfile(input_file):
    print("File does not exist")
    sys.exit(1)

if len(sys.argv) == 3:
    output_file = sys.argv[2]

gbk = SeqIO.to_dict(SeqIO.parse(input_file, "genbank"))
SeqIO.write(gbk.values(), output_file, "fasta")
