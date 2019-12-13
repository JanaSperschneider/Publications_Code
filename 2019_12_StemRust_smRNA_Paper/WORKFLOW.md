### Workflow for analyzing the data

##### Trim small RNA adapters using cutadapt 1.16
```
for f in ${datapath}*.fastq.gz
do
    b=$(basename $f)
    name=$(echo $b | cut -f 1 -d '.')
    ext="_trimmed.fastq"
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 18 -M 28 -q 30 --trim-n --discard-untrimmed -o ${datapath}cutadapt/$name$ext $f
done
```

##### Map against rRNA/tRNA contaminants from RFAM
```
for f in ${datapath}cutadapt/*.fastq
do
    b=$(basename $f)
    name=$(echo $b | cut -f 1 -d '.')
    ext="_unaligned_RFAM_RNACentral.fastq"
    log="_unaligned_RFAM_RNACentral.txt"
    echo $f
    bowtie --sam --un $name$ext ${datapath}RFAM/CONTAMINANTS $f temp.sam >> $name$log 2>&1
done
```

##### Run ShortStack 3.8.5 for rust and wheat
```
ShortStack --bowtie_m 100 --bowtie_cores 4 --readfile G1_trimmed_unaligned_RFAM_RNACentral.fastq \
G2_trimmed_unaligned_RFAM_RNACentral.fastq \
G3_trimmed_unaligned_RFAM_RNACentral.fastq \
U1_trimmed_unaligned_RFAM_RNACentral.fastq \
W0U1_trimmed_unaligned_RFAM_RNACentral.fastq \
W0U2_trimmed_unaligned_RFAM_RNACentral.fastq \
W0U3_trimmed_unaligned_RFAM_RNACentral.fastq \
W3I1_trimmed_unaligned_RFAM_RNACentral.fastq \
W3I2_trimmed_unaligned_RFAM_RNACentral.fastq \
W3I3_trimmed_unaligned_RFAM_RNACentral.fastq \
W5I1_trimmed_unaligned_RFAM_RNACentral.fastq \
W5I2_trimmed_unaligned_RFAM_RNACentral.fastq \
W5I3_trimmed_unaligned_RFAM_RNACentral.fastq \
W7I1_trimmed_unaligned_RFAM_RNACentral.fastq \
W7I2_trimmed_unaligned_RFAM_RNACentral.fastq \
W7I3_trimmed_unaligned_RFAM_RNACentral.fastq \
--genomefile WheatGenome/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta --outdir ${outpath}Wheat_smRNA_PredictionsShortStack

ShortStack --bowtie_m 100 --bowtie_cores 4 --readfile G1_trimmed_unaligned_RFAM_RNACentral.fastq \
G2_trimmed_unaligned_RFAM_RNACentral.fastq \
G3_trimmed_unaligned_RFAM_RNACentral.fastq \
U1_trimmed_unaligned_RFAM_RNACentral.fastq \
W0U1_trimmed_unaligned_RFAM_RNACentral.fastq \
W0U2_trimmed_unaligned_RFAM_RNACentral.fastq \
W0U3_trimmed_unaligned_RFAM_RNACentral.fastq \
W3I1_trimmed_unaligned_RFAM_RNACentral.fastq \
W3I2_trimmed_unaligned_RFAM_RNACentral.fastq \
W3I3_trimmed_unaligned_RFAM_RNACentral.fastq \
W5I1_trimmed_unaligned_RFAM_RNACentral.fastq \
W5I2_trimmed_unaligned_RFAM_RNACentral.fastq \
W5I3_trimmed_unaligned_RFAM_RNACentral.fastq \
W7I1_trimmed_unaligned_RFAM_RNACentral.fastq \
W7I2_trimmed_unaligned_RFAM_RNACentral.fastq \
W7I3_trimmed_unaligned_RFAM_RNACentral.fastq \
--genomefile chr_A_B_unassigned.fasta --outdir ${outpath}Rust_smRNA_PredictionsShortStack
```

##### Prediction of repeats in rust using RepeatModeler 1.0.11 
```
BuildDatabase -name chr_A  chr_A.fasta
RepeatModeler -pa 8 -database chr_A
BuildDatabase -name chr_B chr_B.fasta
RepeatModeler -pa 8 -database chr_B
```
##### Remove TEs from proteome with Blast+ 2.9.0 and following the procedure in https://blaxter-lab-documentation.readthedocs.io/en/latest/filter-repeatmodeler-library.html. This resulted in filtered RepeatModeler files RepeatModeler_ChrsA_consensi.fa.classified and RepeatModeler_ChrsB_consensi.fa.classified

##### RepeatMasking with Repeatmasker 4.0.6 
```
cp /apps/repeatmasker/4.0.6/Libraries/RepeatMaskerLib.embl ${outpath}RepeatMasker_RepeatModeler_ChrsA
/apps/repeatmasker/4.0.6/util/buildRMLibFromEMBL.pl ${outpath}RepeatMasker_RepeatModeler_ChrsA/RepeatMaskerLib.embl > ${outpath}RepeatMasker_RepeatModeler_ChrsA/RepeatMaskerLib.fasta

cat RepeatModeler_ChrsA_consensi.fa.classified ${outpath}RepeatMasker_RepeatModeler_ChrsA/RepeatMaskerLib.fasta > ${outpath}RepeatMasker_RepeatModeler_ChrsA/PGT_Repeats.fasta
RepeatMasker -s -q -lib ${outpath}RepeatMasker_RepeatModeler_ChrsA/PGT_Repeats.fasta -dir ${outpath}RepeatMasker_RepeatModeler_ChrsA -xsmall -pa 8 /datastore/spe12g/Pgt_21_0_Data/chromosomes/chr_A.fasta

cat RepeatModeler_ChrsB_consensi.fa.classified repeatmasker.eukaryotes.fa > ${outpath}RepeatMasker_RepeatModeler_ChrsB/PGT_Repeats.fasta
RepeatMasker -s -q -lib ${outpath}RepeatMasker_RepeatModeler_ChrsB/PGT_Repeats.fasta -dir ${outpath}RepeatMasker_RepeatModeler_ChrsB -xsmall -pa 8 /datastore/spe12g/Pgt_21_0_Data/chromosomes/chr_B.fasta
```
