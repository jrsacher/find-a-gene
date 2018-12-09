# -*- coding: utf-8 -*-
"""
Created on Sun 09 2018
Last update Sun 09 2018

@author: jrsacher

Find-a-gene v 0.3
See chapter 4 in Bioinformatics and Functional Genomics (3rd ed)
http://bioinfbook.org/php/C4E3k
"""

import csv
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from distutils.util import strtobool

# Change this if you aren't me!
Entrez.email = "joshuasacher@gmail.com"

# Open file, only keep novel genes
with open("novel_genes.csv", "r") as infile:
    reader = csv.DictReader(infile)
    novels = [row for row in reader if strtobool(row["novel"].lower())]

gis = ["gi|" + gene["gi"].replace(" ***", "") for gene in novels]

with open("html.txt", "w") as f:
    f.write('<table class="table table-hover">\n')
    f.write('\t<thead class="thead-dark">\n')
    f.write("\t\t<th>ID</th>\n")
    f.write("\t\t<th>Source</th>\n")
    f.write("\t\t<th>Description</th>\n")
    f.write("\t\t<th>Link</th>\n")
    f.write("\t</thead>\n")
    f.write("\t<tbody>\n")

    with Entrez.efetch(db="nucest", rettype="gb", retmode="text", id=gis) as handle:
        seqs = list(SeqIO.parse(handle, "gb"))

    link = "https://www.ncbi.nlm.nih.gov/nucest/"
    for seq in seqs:
        f.write("\t\t<tr>\n")
        f.write(f'\t\t\t<td>{seq.id}</td>\n')
        f.write(f'\t\t\t<td>{seq.annotations["source"]}</td>\n')
        f.write(f'\t\t\t<td>{seq.description}</td>\n')
        f.write(f'\t\t\t<td><a href="{link + seq.id}">info</a></td>\n')
        f.write("\t\t</tr>\n")
    f.write("\t</tbody>\n")
    f.write("</table>\n")
