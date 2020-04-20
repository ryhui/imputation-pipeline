"""
21th January 2019

pmd_filter.py

The script filters vcf file to remove sites where calls could be down to deamination. 
Adapted from Lara Cassidy's filter_vcf_for_imputation_v3.py

Running the script:
python pmd_filter.py input.vcf[.gz] | bgzip -c > output.vcf.gz
"""

import sys
import os
import gzip

args = sys.argv
if args[1].endswith('.gz'):
    input_file = gzip.open(args[1], 'r')
else:
    input_file = open(args[1], 'rb')

for each in input_file:
    each = each.decode('utf-8').strip()
    col = each.split('\t')

    if each.startswith('#'):
        print(each)
    else:
        info = col[9].split(":")

        format = col[8].split(':')
        genotype = info[format.index('GT')]
        likelihood = info[format.index('PL')]

        pl = [int(p) for p in likelihood.split(',')]

        if ((col[3] == "G" and col[4] == "A") or (col[3] == "C" and col[4] == "T")) and (pl[2] == 0 or pl[1] == 0): 
            continue
        elif ((col[3] == "A" and col[4] == "G") or (col[3] == "T" and col[4] == "C")) and (pl[0] == 0 or pl[1] == 0):
            continue
        else:
            print(each)

input_file.close()
