import os 
import sys
import operator
# -----------------------------------------------------------------------------------------------------------
def collect_hits_from_mapped_reads(SAMFILE):
    hits = {}

    with open(SAMFILE) as f: 
        content = f.readlines()
        for line in content:    
            # This is where the read mapping starts
            read = line.split('\t')[0]
            hit_position = int(line.split('\t')[3])
            chromosome = line.split('\t')[2]
            if chromosome in hits:
                previous_positions = hits[chromosome]
                hits[chromosome] = previous_positions + [(hit_position, read)]
            else:
                hits[chromosome] = [(hit_position, read)]

    return hits
# -----------------------------------------------------------------------------------------------------------
def sRNA_loci_read_coverage(hits, identifiers, total_reads):

    best_hits = {}
    
    for sRNA in identifiers:
        chromosome = sRNA.split('locus:')[1].split(':')[0]
        locus = sRNA.split('locus:')[1].split(':')[1].split('|')[0]
        start = int(locus.split('-')[0])
        end = int(locus.split('-')[1])
        
        sRNA_read_coverage = 0

        if chromosome in hits:

            read_positions = hits[chromosome]
            reads_covered = []

            for pos, read in read_positions:
                if pos >= start and pos <= end:
                    sRNA_read_coverage += 1
                    reads_covered.append(read)

        if sRNA_read_coverage > 0:
            best_hits[sRNA] = len(list(set(reads_covered)))#sRNA_read_coverage

    return best_hits
# -----------------------------------------------------------------------------------------------------------
HAPLOTYPE_MAPPING = {'chr1A': ['chr1B'],
                     'chr2A': ['chr2B'],
                     'chr3A': ['chr3B', 'chr5B'],
                     'chr4A': ['chr4B'],
                     'chr5A': ['chr5B', 'chr3B'],
                     'chr6A': ['chr6B'],
                     'chr7A': ['chr7B'],
                     'chr8A': ['chr8B', 'chr16B'],
                     'chr9A': ['chr9B'],
                     'chr10A': ['chr10B'],
                     'chr11A': ['chr11B'],
                     'chr12A': ['chr12B'],
                     'chr13A': ['chr13B'],
                     'chr14A': ['chr14B'],
                     'chr15A': ['chr15B'],
                     'chr16A': ['chr16B', 'chr8B'],
                     'chr17A': ['chr17B'],
                     'chr18A': ['chr18B'],
                     'chr1B': ['chr1A'],
                     'chr2B': ['chr2A'],
                     'chr3B': ['chr3A', 'chr5A'],
                     'chr4B': ['chr4A'],
                     'chr5B': ['chr5A', 'chr3A'],
                     'chr6B': ['chr6A'],
                     'chr7B': ['chr7A'],
                     'chr8B': ['chr8A', 'chr16A'],
                     'chr9B': ['chr9A'],
                     'chr10B': ['chr10A'],
                     'chr11B': ['chr11A'],
                     'chr12B': ['chr12A'],
                     'chr13B': ['chr13A'],
                     'chr14B': ['chr14A'],
                     'chr15B': ['chr15A'],
                     'chr16B': ['chr16A', 'chr8A'],
                     'chr17B': ['chr17A'],
                     'chr18B': ['chr18A']
                     }                     
# -----------------------------------------------------------------------------------------------------------
FASTA_FILE = sys.argv[1]
RESULTS_FILE = '../' + FASTA_FILE.replace('.fasta', '_homologous_smRNAs.csv')
o = open(RESULTS_FILE, 'w')

# Get the smRNA sequences and identifiers
identifiers = []
sequences = []

with open('/scratch1/spe12g/StemRust_smRNA/Rust_smRNA_PredictionsShortStack/' + FASTA_FILE) as f: 
    content = f.readlines()

    for position, line in enumerate(content):
        if '>' in line:
            identifiers.append(line)
        else:
            sequences.append(line.strip())
# -----------------------------------------------------------------------------------------------------------
identifiers_chrs = {}
for ident, seq in zip(identifiers, sequences):
    chromosome = ident.split('locus:')[1].split(':')[0]
    identifiers_chrs[ident] = seq

print 'All small RNAs:', len(identifiers_chrs)
count_A = sum([1 for ident in identifiers if 'A' in ident.split('locus:')[1].split(':')[0]])
count_B = sum([1 for ident in identifiers if 'B' in ident.split('locus:')[1].split(':')[0]])
print 'A chromosomes:', count_A, 100.0*count_A/float(count_A + count_B)
print 'B chromosomes:', count_B, 100.0*count_B/float(count_A + count_B)

# -----------------------------------------------------------------------------------------------------------
count_sRNAs_with_homolog, count_sRNAs_without_homolog = 0, 0
homolog_on_other_haplotype, homolog_on_same_haplotype = 0, 0
sRNAs_on_chromosomes = 0

for smRNA_chr, seq in identifiers_chrs.items():
    chromosome = smRNA_chr.split('locus:')[1].split(':')[0]
    locus = smRNA_chr.split('locus:')[1].split(':')[1].split('|')[0]
    start = locus.split('-')[0]
    end = locus.split('-')[1]

    if 'tig' in chromosome:
        pass
    else:
        sRNAs_on_chromosomes += 1
        print '---------------------'
        print smRNA_chr
        print start, end
        print
        # Map all the reads from that region to all chromosomes to find the homologous region
        os.popen('samtools view /scratch1/spe12g/StemRust_smRNA/Rust_smRNA_PredictionsShortStack/merged_alignments.bam ' + chromosome + ':' + start + '-' + end + ' > test.sam')
        os.popen(''' awk '/^@/{ next; } { print ">"$1; print $10 }' test.sam > test.fasta''')
        os.popen('bowtie -f -v2 -a -m 100 --best --strata --sam /scratch1/spe12g/Bowtie_Indices/chr_A_B_unassigned test.fasta > test.sam')
        os.popen('samtools view -b -o test.bam test.sam')
        os.popen('bedtools coverage -nonamecheck -a Results.bed -b test.bam > out.txt')
        # -----------------------------------------------------------------------------------------------------------
        does_sRNA_have_homolog = False
        candidates = []

        # Go through all the sRNA loci and look at their read coverage 
        with open('out.txt') as f: 
            content = f.readlines()        
            for line in content:
                possible_homolog_chromosome, possible_homolog_start, possible_homolog_end = line.split('\t')[0], line.split('\t')[1], line.split('\t')[2]
                read_coverage = float(line.split('\t')[6])
                # Do not look at the smRNA origin locus
                if possible_homolog_chromosome == chromosome and possible_homolog_start == start and possible_homolog_end == end:
                    pass
                else:
                    # ----------------------------------------------
                    if 100.0*read_coverage > 25.0: # The fraction of bases in the sRNA that had non-zero coverage
                        # ----------------------------------------------
                        candidates.append((100.0*read_coverage, possible_homolog_chromosome, possible_homolog_start, possible_homolog_end))

        if candidates:
            best_homologous_sRNA = sorted(candidates)[-1]
            homolog_chromosome, homolog_start, homolog_end = best_homologous_sRNA[1], best_homologous_sRNA[2], best_homologous_sRNA[3]
            read_coverage = best_homologous_sRNA[0]

            print
            print 'Best possible homolog', homolog_chromosome, homolog_start, homolog_end, 'with read coverage', 100.0*read_coverage
            print 
            # ----------------------------------------------
            output = smRNA_chr.strip() + ',' + chromosome + ',' + start + ',' + end + ',' + homolog_chromosome + ',' + str(homolog_start) + ',' + str(homolog_end)
            o.writelines(output + '\n')
            # ----------------------------------------------
            other_haplotype = HAPLOTYPE_MAPPING[chromosome]

            if homolog_chromosome in other_haplotype:
                homolog_on_same_haplotype += 1
            else:
                homolog_on_other_haplotype += 1
            
            does_sRNA_have_homolog = True
                    
        if does_sRNA_have_homolog == False:
            count_sRNAs_without_homolog += 1
        else:
            count_sRNAs_with_homolog += 1
        print '---------------------'
# -----------------------------------------------------------------------------------------------------------
o.close()
# -----------------------------------------------------------------------------------------------------------
print RESULTS_FILE
print 
print 'All small RNAs on chromosomes:', sRNAs_on_chromosomes
print 'sRNAs with homologous counterpart', count_sRNAs_with_homolog, 100.0*count_sRNAs_with_homolog/float(sRNAs_on_chromosomes)
print 'sRNAs without homologous counterpart', count_sRNAs_without_homolog, 100.0*count_sRNAs_without_homolog/float(sRNAs_on_chromosomes)
print
print 'Homologous counterpart on corresponding haplotype', homolog_on_same_haplotype, 100.0*homolog_on_same_haplotype/float(homolog_on_same_haplotype + homolog_on_other_haplotype)
print 'Homologous counterpart on different haplotype', homolog_on_other_haplotype, 100.0*homolog_on_other_haplotype/float(homolog_on_same_haplotype + homolog_on_other_haplotype)
print