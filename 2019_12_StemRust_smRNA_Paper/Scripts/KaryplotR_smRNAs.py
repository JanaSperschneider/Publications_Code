import os
import sys
# --------------------------------------------------------------
fileinput = sys.argv[1]
f = open(fileinput, 'r')
content = f.readlines()
f.close()

output = open(fileinput.replace('.fasta','_karyplotR.txt'), 'w')
# --------------------------------------------------------------
# --------------------------------------------------------------
# --------------------------------------------------------------
for line in content:
    if '>' in line:
        cluster = line.split('|')[-1].strip()        
        chromosome = line.split('|locus:')[1].split(':')[0]
        start = line.split('|locus:')[1].split(':')[1].split('-')[0]
        end = line.split('|locus:')[1].split(':')[1].split('-')[1].split('|')[0]
        output.writelines(chromosome + '\t' + start + '\t' + end + '\n')
output.close()