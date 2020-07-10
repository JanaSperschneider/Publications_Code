import os
import sys
# --------------------------------------------------------------
fileinput = sys.argv[1]
f = open(fileinput, 'r')
content = f.readlines()
f.close()

align_file = fileinput.replace('Results.fasta', 'merged_alignments.bam')
# --------------------------------------------------------------
total_count, count_A, count_C, count_G, count_U = 0, 0, 0, 0, 0
seq_list = []
# --------------------------------------------------------------
# --------------------------------------------------------------
for line in content:
    if '>' in line:
        chromosome = line.split('|locus:')[1].split(':')[0]
        start = line.split('|locus:')[1].split(':')[1].split('-')[0]
        end = line.split('|locus:')[1].split(':')[1].split('-')[1].split('|')[0]
        # Now go through SAM file and collect all the other reads
        os.popen('samtools view ' + align_file + ' ' + chromosome + ':' + start + '-' + end + ' > test.sam')


        reads = []

        f = open('test.sam', 'r')
        samfile = f.readlines()
        f.close()

        for entry in samfile:
            reads.append(entry.split('\t')[9])

        for read in reads:
            if read[0] == 'A':
                count_A += 1
            if read[0] == 'C':
                count_C += 1
            if read[0] == 'G':
                count_G += 1
            if read[0] == 'T':
                count_U += 1
            total_count += 1
        seq_list.append(read.strip())

print total_count, 'sequences analyzed.'
print count_A, count_C, count_G, count_U, count_A+count_C+count_G+count_U == total_count
print 'A', 100.0*count_A/float(total_count), 
print 'C', 100.0*count_C/float(total_count), 
print 'G', 100.0*count_G/float(total_count), 
print 'U', 100.0*count_U/float(total_count)
print 100.0*count_A/float(total_count) + 100.0*count_C/float(total_count) + 100.0*count_G/float(total_count) + 100.0*count_U/float(total_count) == 100.0
print

print "% of reads with length, %U, %A, %C, %G" 
for i in xrange(18, 29):
    print str(i), round(100.0*sum([1 for seq in seq_list if len(seq) == float(i)])/len(seq_list),1),  
    if 100.0*sum([1 for seq in seq_list if len(seq) == float(i)])/len(seq_list) != 0.0:
        print ' %A:', round(100.0*sum([1 for seq in seq_list if seq[0] == 'A' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
        print ' %C:', round(100.0*sum([1 for seq in seq_list if seq[0] == 'C' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
        print ' %G:', round(100.0*sum([1 for seq in seq_list if seq[0] == 'G' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
        print ' %U:', round(100.0*sum([1 for seq in seq_list if seq[0] == 'T' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1)
    else:
        print ' %A: 0.0', 
        print ' %C: 0.0', 
        print ' %G: 0.0', 
        print ' %U: 0.0'

print "Number of reads with length"
for i in xrange(18, 29):
    print str(i), sum([1 for seq in seq_list if len(seq) == float(i)])

#------------------------------------------------------
print "% of reads with length"
for i in xrange(18, 29):
    print round(100.0*sum([1 for seq in seq_list if len(seq) == float(i)])/len(seq_list),1), 
    print ',', 
#------------------------------------------------------
print
i=21
print 'L21 <- c(',
if 100.0*sum([1 for seq in seq_list if len(seq) == float(i)])/len(seq_list) != 0.0:
    print round(100.0*sum([1 for seq in seq_list if seq[0] == 'A' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
    print round(100.0*sum([1 for seq in seq_list if seq[0] == 'C' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
    print round(100.0*sum([1 for seq in seq_list if seq[0] == 'G' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
    print round(100.0*sum([1 for seq in seq_list if seq[0] == 'T' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1)
else:
    print '0.0', 

#------------------------------------------------------
print
i=22
print 'L22 <- c(',
if 100.0*sum([1 for seq in seq_list if len(seq) == float(i)])/len(seq_list) != 0.0:
    print round(100.0*sum([1 for seq in seq_list if seq[0] == 'A' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
    print round(100.0*sum([1 for seq in seq_list if seq[0] == 'C' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
    print round(100.0*sum([1 for seq in seq_list if seq[0] == 'G' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1),
    print round(100.0*sum([1 for seq in seq_list if seq[0] == 'T' and len(seq) == float(i)])/float(sum([1 for seq in seq_list if len(seq) == float(i)])),1)
else:
    print '0.0', 
