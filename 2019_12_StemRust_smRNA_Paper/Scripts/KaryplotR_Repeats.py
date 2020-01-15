import os
import sys
import subprocess
import operator
import re
# -----------------------------------------------------------------------------------------------------------
'''
Load the repeats into a dictionary
'''
f = open('../chr_A_B.fasta.out', 'r')
repeatcontent = f.readlines()
f.close()
# -----------------------------------------------------------------------------------------------------------
f = open('../chrs_repeatmasker_no_simple_repeats.txt', 'w')

# Load repeat dictionary
for line in repeatcontent[3:]:
    line = line.split(' ')
    line = [entry for entry in line if entry != '']
    contig = line[4]
    contig = contig.replace('_','')

    start = int(line[5])
    end = int(line[6])
    length = end-start

    repeat = line[9]
    repeat_class = line[10]

    if repeat_class != 'Simple_repeat' and repeat_class != 'Low_complexity' and not 'tig' in contig:

        output = contig + '\t' + str(start) + '\t' + str(end) + '\n'
        f.writelines(output)

f.close()
# -----------------------------------------------------------------------------------------------------------
f = open('../chrs_repeatmasker_gypsys.txt', 'w')

# Load repeat dictionary
for line in repeatcontent[3:]:
    line = line.split(' ')
    line = [entry for entry in line if entry != '']
    contig = line[4]
    contig = contig.replace('_','')

    start = int(line[5])
    end = int(line[6])
    length = end-start

    repeat = line[9]
    repeat_class = line[10]  

    if 'LTR/Gypsy' in repeat_class and not 'tig' in contig:

        output = contig + '\t' + str(start) + '\t' + str(end) + '\n'
        f.writelines(output)

f.close()
# -----------------------------------------------------------------------------------------------------------
f = open('../chrs_repeatmasker_DNAtransposons.txt', 'w')

# Load repeat dictionary
for line in repeatcontent[3:]:
    line = line.split(' ')
    line = [entry for entry in line if entry != '']
    contig = line[4]
    contig = contig.replace('_','')
    start = int(line[5])
    end = int(line[6])
    length = end-start

    repeat = line[9]
    repeat_class = line[10].split('/')[0]     

    if 'DNA' in repeat_class and not 'tig' in contig:

        output = contig + '\t' + str(start) + '\t' + str(end) + '\n'
        f.writelines(output)

f.close()