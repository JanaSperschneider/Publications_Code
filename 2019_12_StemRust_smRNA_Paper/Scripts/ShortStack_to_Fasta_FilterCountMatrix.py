import os
import sys
###############################################################
RPM_THRESHOLD = 5.0
###############################################################
fileinput = sys.argv[1]

f = open(fileinput, 'r')
content = f.readlines()
f.close()
###############################################################
f = open(fileinput.replace('.txt','.fasta'), 'w')

# For each line, extract the small RNA, there will be duplicate entries, so preserve the genomic location
small_RNAs, identifiers, phase_scores, rpms, clusters = [], [], [], [], []
count_invalid = 0

#Locus	Name	Length	Reads	RPM	UniqueReads	FracTop	Strand	MajorRNA	MajorRNAReads	Complexity	DicerCall	MIRNA	PhaseScore	Short	Long	20	21	22	23	24

for line in content[1:]:

    locus = line.split('\t')[0]
    cluster = line.split('\t')[1]
    rpm = float(line.split('\t')[4])
    smallRNA = line.split('\t')[8]
    dicercall = line.split('\t')[11]
    majorRNAreads = int(line.split('\t')[9])
    phasing = line.split('\t')[13]

    if dicercall != 'N' and rpm >= RPM_THRESHOLD:

        small_RNAs.append(smallRNA)
        identifiers.append(locus)
        phase_scores.append(phasing)
        rpms.append(rpm)
        clusters.append(cluster)
    else:
        count_invalid += 1

for index, smallRNA in enumerate(small_RNAs):

    f.writelines('>smRNA_' + str(index) + '|rpm:' + str(rpms[index]) + '|locus:' + identifiers[index] + '|' + clusters[index] + '\n')
    f.writelines(str(smallRNA.replace('U','T')) + '\n')

f.close()
###############################################################
print fileinput
print "Number of total small RNA sequences:", len(small_RNAs)
print "Number of unique small RNA sequences:", len(list(set(small_RNAs)))
print "Removed because of invalid dicer calls:", count_invalid
###############################################################
f = open(fileinput.replace('Results.txt','Counts.txt'), 'r')
content = f.readlines()
f.close()
###############################################################
f = open(fileinput.replace('Results.txt','Counts_Filtered.txt'), 'w')
f.writelines(content[0])

count = 0
for line in content[1:]:
    cluster = line.split('\t')[1]
    for smallRNA in clusters:        
        cluster_smRNA = smallRNA.split('|')[-1].strip()
        if cluster == cluster_smRNA:
            f.writelines(line)
            count += 1   
            break
f.close()

print count
