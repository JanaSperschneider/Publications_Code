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

##### Genomic origins of sRNAs with bedtools 2.28.0
```
# Turn the gene annotation file into bed file format
awk -F'\t' -v OFS='\t' 'NR>=2{ if ($3 == "gene") { print } }' Puccinia_graminis_tritici_21-0.gff3 | awk '{print $1 "\t" $4 "\t" $5 "\t" $9 "\t" 0 "\t" $7}' > Puccinia_graminis_tritici_21-0.bed
sortBed -i Puccinia_graminis_tritici_21-0.bed > Puccinia_graminis_tritici_21-0.sorted.bed

# Turn the repeat annotation file from RepeatMasker into bed file format
cat chr_A_B.fasta.out | awk '{print $5 "\t" $6 "\t" $7 "\t" $10 "\t" $11}' > chr_A_B.fasta.out.bed
awk -F'\t' -v OFS='\t' 'NR>=3{sub(/_/, "", $1)} 1' chr_A_B.fasta.out.bed | tail -n+4 > f1.txt
sortBed -i f1.txt > chr_A_B.fasta.out.sorted.bed
rm f1.txt

# Turn the sRNAs locations into bed file format
grep '>' ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.fasta > smRNAs_up_late_infection.IDs
cat smRNAs_up_late_infection.IDs | cut -d ':' -f3,4 | cut -d '|' -f1 | sed 's/:/\t/g' | sed 's/-/\t/g' > smRNAs_up_late_infection.bed
sortBed -i smRNAs_up_late_infection.bed > smRNAs_up_late_infection.sorted.bed
grep '>' ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_upSpores.fasta > smRNAs_upSpores.IDs
cat smRNAs_upSpores.IDs | cut -d ':' -f3,4 | cut -d '|' -f1 | sed 's/:/\t/g' | sed 's/-/\t/g' > smRNAs_upSpores.bed
sortBed -i smRNAs_upSpores.bed > smRNAs_upSpores.sorted.bed
grep '>' ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_noDE.fasta > smRNAs_noDE.IDs
cat smRNAs_noDE.IDs | cut -d ':' -f3,4 | cut -d '|' -f1 | sed 's/:/\t/g' | sed 's/-/\t/g' > smRNAs_noDE.bed
sortBed -i smRNAs_noDE.bed > smRNAs_noDE.sorted.bed
grep '>' ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_early_infection.fasta > smRNAs_up_early_infection.IDs
cat smRNAs_up_early_infection.IDs | cut -d ':' -f3,4 | cut -d '|' -f1 | sed 's/:/\t/g' | sed 's/-/\t/g' > smRNAs_up_early_infection.bed
sortBed -i smRNAs_up_early_infection.bed > smRNAs_up_early_infection.sorted.bed

bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_late_infection.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo > smRNAs_up_late_infection_repeats_overlapping.bed
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_late_infection.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo | awk '{print $1 "\t" $2 "\t" $3}' | uniq | wc -l
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_late_infection.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo | awk '{print $8}' | sort | uniq -c | sort -g

bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_upSpores.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo > smRNAs_upSpores_repeats_overlapping.bed
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_upSpores.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo | awk '{print $1 "\t" $2 "\t" $3}' | uniq | wc -l
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_upSpores.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo | awk '{print $8}' | sort | uniq -c | sort -g

bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_early_infection.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo > smRNAs_up_early_infection_repeats_overlapping.bed
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_early_infection.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo | awk '{print $1 "\t" $2 "\t" $3}' | uniq | wc -l
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_early_infection.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo | awk '{print $8}' | sort | uniq -c | sort -g

bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_noDE.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo > smRNAs_noDE_repeats_overlapping.bed
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_noDE.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo | awk '{print $1 "\t" $2 "\t" $3}' | uniq | wc -l
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_noDE.sorted.bed -b chr_A_B.fasta.out.sorted.bed -wo | awk '{print $8}' | sort | uniq -c | sort -g

# sRNAs overlapping/containing genes 
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_upSpores.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed -wo > smRNAs_upSpores_genes_overlapping.bed
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_upSpores.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed -wo | awk '{print $1 "\t" $2 "\t" $3}' | uniq | wc -l
cat smRNAs_upSpores_genes_overlapping.bed | awk '{print $7}' | uniq | sed 's/ID=//g' | cut -d ';' -f1 > smRNAs_upSpores_genes_overlapping.gene.IDs

bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_late_infection.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed -wo > smRNAs_up_late_infection_genes_overlapping.bed
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_late_infection.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed -wo | awk '{print $1 "\t" $2 "\t" $3}' | uniq | wc -l
cat smRNAs_up_late_infection_genes_overlapping.bed | awk '{print $7}' | uniq | sed 's/ID=//g' | cut -d ';' -f1 > smRNAs_up_late_infection_genes_overlapping.gene.IDs

bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_early_infection.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed -wo > smRNAs_up_early_infection_genes_overlapping.bed
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_up_early_infection.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed -wo | awk '{print $1 "\t" $2 "\t" $3}' | uniq | wc -l
cat smRNAs_up_early_infection_genes_overlapping.bed | awk '{print $7}' | uniq | sed 's/ID=//g' | cut -d ';' -f1 > smRNAs_up_early_infection_genes_overlapping.gene.IDs

bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_noDE.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed -wo > smRNAs_noDE_genes_overlapping.bed
bedtools intersect -nonamecheck -f 0.25 -F 0.25 -a smRNAs_noDE.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed -wo | awk '{print $1 "\t" $2 "\t" $3}' | uniq | wc -l
cat smRNAs_noDE_genes_overlapping.bed | awk '{print $7}' | uniq | sed 's/ID=//g' | cut -d ';' -f1 > smRNAs_noDE_genes_overlapping.gene.IDs
```

##### Figure 9: TEs targeted by sRNAs are associated with reduced expression of overlapping genes using bowtie 1.1.2
```
bowtie -f -v0 -a --best --strata --sam chr_A_B_unassigned ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.fasta ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.sam
sam2bed < ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.sam > ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.bed
# Bed file of the sRNA mappings back to the genome
sortBed -i ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.bed > ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.sorted.bed
cat ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.sorted.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 0 "\t" $6}' > ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.bed

# Get the repeats that do and don't overlap with sRNAs
bedtools intersect -nonamecheck -a chr_A_B.fasta.out.sorted.bed -b ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.bed -wo > ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection_repeats_overlapping.bed
bedtools intersect -nonamecheck -v -a chr_A_B.fasta.out.sorted.bed -b ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.bed -wo > ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection_repeats_not_overlapping.bed

# Now get the repeats 
cat ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection_repeats_overlapping.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > ${outpath}/Rust_smRNA_PredictionsShortStack/repeats_overlapping_smRNAs_up_late_infection.bed
sortBed -i ${outpath}/Rust_smRNA_PredictionsShortStack/repeats_overlapping_smRNAs_up_late_infection.bed > ${outpath}/Rust_smRNA_PredictionsShortStack/repeats_overlapping_smRNAs_up_late_infection.sorted.bed
cat ${outpath}/Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection_repeats_not_overlapping.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > ${outpath}/Rust_smRNA_PredictionsShortStack/repeats_not_overlapping_smRNAs_up_late_infection.bed
sortBed -i ${outpath}/Rust_smRNA_PredictionsShortStack/repeats_not_overlapping_smRNAs_up_late_infection.bed > ${outpath}/Rust_smRNA_PredictionsShortStack/repeats_not_overlapping_smRNAs_up_late_infection.sorted.bed

# Now get the genes that overlap with those repeats
bedtools closest -nonamecheck -d -a ${outpath}/Rust_smRNA_PredictionsShortStack/repeats_overlapping_smRNAs_up_late_infection.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed > repeats_overlapping_smRNAs_up_late_infection_closest_genes.bed
bedtools closest -nonamecheck -d -a ${outpath}/Rust_smRNA_PredictionsShortStack/repeats_not_overlapping_smRNAs_up_late_infection.sorted.bed -b Puccinia_graminis_tritici_21-0.sorted.bed > repeats_not_overlapping_smRNAs_up_late_infection_closest_genes.bed
```
