import os
import sys
###############################################################
RPM_THRESHOLD = 5.0
###############################################################
fileinput = sys.argv[1]
###############################################################
f = open(fileinput, 'r')
content = f.readlines()
f.close()
###############################################################
# For each line, extract the small RNA, there will be duplicate entries, so preserve the genomic location
small_RNAs = {}
siRNAs, miRNAs = {}, {}
count_invalid = 0
###############################################################
#Locus	Name	Length	Reads	RPM	UniqueReads	FracTop	Strand	MajorRNA	MajorRNAReads	Complexity	DicerCall	MIRNA	PhaseScore	Short	Long	20	21	22	23	24
###############################################################
f = open(fileinput.replace('.txt','.fasta'), 'w')

f_siRNA = open(fileinput.replace('.txt','_siRNAs.fasta'), 'w')
f_miRNA = open(fileinput.replace('.txt','_miRNAs.fasta'), 'w')

for line in content[1:]:
    locus = line.split('\t')[0]
    cluster = line.split('\t')[1]
    rpm = float(line.split('\t')[4])
    sequence = line.split('\t')[8]
    majorRNAreads = int(line.split('\t')[9])    
    dicercall = line.split('\t')[11]
    miRNA = line.split('\t')[12]
    phasing = line.split('\t')[13]

    if dicercall != 'N' and rpm >= RPM_THRESHOLD:        
        small_RNAs[cluster] = [rpm, locus, sequence]

        if miRNA == 'Y':
            miRNAs[cluster] = [rpm, locus, sequence]
        else:
            siRNAs[cluster] = [rpm, locus, sequence]

    else:
        count_invalid += 1

for cluster, values in small_RNAs.items():
    rpm, locus, sequence = values[0], values[1], values[2]
    f.writelines('>smRNA' + '|rpm:' + str(rpm) + '|locus:' + locus + '|' + cluster + '\n')
    f.writelines(str(sequence.replace('U','T')) + '\n')
f.close()

for cluster, values in siRNAs.items():
    rpm, locus, sequence = values[0], values[1], values[2]
    f_siRNA.writelines('>smRNA' + '|rpm:' + str(rpm) + '|locus:' + locus + '|' + cluster + '\n')
    f_siRNA.writelines(str(sequence.replace('U','T')) + '\n')
f_siRNA.close()

for cluster, values in miRNAs.items():
    rpm, locus, sequence = values[0], values[1], values[2]
    f_miRNA.writelines('>smRNA' + '|rpm:' + str(rpm) + '|locus:' + locus + '|' + cluster + '\n')
    f_miRNA.writelines(str(sequence.replace('U','T')) + '\n')
f_miRNA.close()

print fileinput
print "Number of total small RNA sequences:", len(small_RNAs)
print "Removed because of invalid dicer calls and other criteria:", count_invalid
###############################################################
print
print "Number of siRNAs: ", len(siRNAs)
print "Number of miRNAs: ", len(miRNAs)
