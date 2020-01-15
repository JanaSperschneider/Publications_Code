#### HiC-Pro was run to produce contact maps and pinpoint the location of the rust centromeres

```
module load hic-pro/2.11.1

# Generate BED file of the restriction fragments after digestion.
# Sau3A [^GATC]
# -------------------------------------------------------------
/apps/hic-pro/2.11.1/HiC-Pro-master/bin/utils/digest_genome.py -r ^GATC -o chr_A_B_unassigned.bed chr_A_B_unassigned.fasta
# -------------------------------------------------------------
module load bowtie/2.2/9
bowtie2-build chr_A_B_unassigned.fasta chr_A_B_unassigned
# -------------------------------------------------------------
cat Pgt_HiCPro/rawdata/PgtHIC_S1_L001_R1_001.fastq.gz Pgt_HiCPro/rawdata/PgtHIC_S1_L002_R1_001.fastq.gz Pgt_HiCPro/rawdata/PgtHIC_S1_L003_R1_001.fastq.gz Pgt_HiCPro/rawdata/PgtHIC_S1_L004_R1_001.fastq.gz > Pgt_HiCPro/rawdata/PgtHIC_R1.fastq.gz
cat Pgt_HiCPro/rawdata/PgtHIC_S1_L001_R2_001.fastq.gz Pgt_HiCPro/rawdata/PgtHIC_S1_L002_R2_001.fastq.gz Pgt_HiCPro/rawdata/PgtHIC_S1_L003_R2_001.fastq.gz Pgt_HiCPro/rawdata/PgtHIC_S1_L004_R2_001.fastq.gz > Pgt_HiCPro/rawdata/PgtHIC_R2.fastq.gz

HiC-Pro -i Pgt_HiCPro/ -o Pgt_HiCPro_Results/ -c config-hicpro.txt
 -------------------------------------------------------------
## Plot the genome-wide map
module load hicexplorer/3.0.1

hicConvertFormat --matrices Pgt_HiCPro_Results/hic_results/matrix/rawdata/iced/20000/PgtHIC_chr_A_B_unassigned.bwt2pairs_20000_iced.matrix --inputFormat hicpro --outputFormat h5 --outFileName 20000_iced.matrix --bedFileHicpro Pgt_HiCPro_Results/hic_results/matrix/rawdata/raw/20000/PgtHIC_chr_A_B_unassigned.bwt2pairs_20000_abs.bed

hicPlotMatrix --matrix 20000_iced.matrix.h5 --out chr1A.png --dpi 600 --log1p --perChromosome --chr chr1A
hicPlotMatrix --matrix 20000_iced.matrix.h5 --out chr1A_Centromere.png --dpi 600 --log1p --region chr1A:2000000-2800000 
```


#### config-hicpro.txt
```
# Please change t	he variable settings below if necessary

#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !!
#######################################################################
N_CPU = 2
LOGFILE = hicpro.log

JOB_NAME = 
JOB_MEM = 
JOB_WALLTIME = 
JOB_QUEUE = 
JOB_MAIL = 

#########################################################################
## Data
#########################################################################

PAIR1_EXT = _R1
PAIR2_EXT = _R2

#######################################################################
## Alignment options
#######################################################################

MIN_MAPQ = 10

BOWTIE2_IDX_PATH = HiCPro
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = chr_A_B_unassigned
GENOME_SIZE = chr_A_B_unassigned.lengths

#######################################################################
## Allele specific analysis
#######################################################################

ALLELE_SPECIFIC_SNP = 

#######################################################################
## Capture Hi-C analysis
#######################################################################

CAPTURE_TARGET =
REPORT_CAPTURE_REPORTER = 1

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = chr_A_B_unassigned.bed
LIGATION_SITE =
MIN_FRAG_SIZE = 
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 20000 40000 150000 500000 1000000
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1

```
