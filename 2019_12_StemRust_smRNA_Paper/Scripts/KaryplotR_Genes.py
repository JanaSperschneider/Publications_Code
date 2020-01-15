import os
import sys
import subprocess
import operator
import re
# -----------------------------------------------------------------------------------------------------------
'''
Load the genes into a dictionary
'''
f = open('../Puccinia_graminis_tritici_21-0.gff3', 'r')
genecontent = f.readlines()
f.close()
# -----------------------------------------------------------------------------------------------------------
f = open('../Genes_KaryoplotR.txt', 'w')

for line in genecontent[1:]:
    line = line.split('\t')
    contig = line[0]
    entry_type = line[2]
    start = int(line[3])
    end = int(line[4])

    if entry_type == 'gene' and not 'tig' in contig:
        output = contig + '\t' + str(start) + '\t' + str(end) + '\n'
        f.writelines(output)

f.close()
# -----------------------------------------------------------------------------------------------------------
f = open('../Tx2_Abundance.txt', 'r')
content = f.readlines()
f.close()

EXPRESSED = {}

for line in content[1:]:
    gene = line.split('\t')[0].replace('"','')
    expression = sum([float(entry) for entry in line.split('\t')[1:]])
    gene = gene.split('-T')[0]
    if gene in EXPRESSED:
        EXPRESSED[gene] += expression
    else:
        EXPRESSED[gene] = expression
# -----------------------------------------------------------------------------------------------------------
f = open('../ExpressedGenes_KaryoplotR.txt', 'w')

for line in genecontent[1:]:
    line = line.split('\t')
    contig = line[0]
    entry_type = line[2]
    start = int(line[3])
    end = int(line[4])
    gene = line[8].split(';')[0]
    gene = gene.replace('ID=','')

    if entry_type == 'gene' and not 'tig' in contig:
        if gene in EXPRESSED:
            if EXPRESSED[gene] > 0.0:
                output = contig + '\t' + str(start) + '\t' + str(end) + '\n'
                f.writelines(output)

f.close()
# -----------------------------------------------------------------------------------------------------------
