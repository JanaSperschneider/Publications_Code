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
ShortStack --bowtie_m all --bowtie_cores 4 --readfile G1_trimmed_unaligned_RFAM_RNACentral.fastq \
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
--genomefile WheatGenome/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta --outdir ${outpath}Wheat_smRNA_Predictions_ShortStack

ShortStack --bowtie_m all --bowtie_cores 4 --readfile G1_trimmed_unaligned_RFAM_RNACentral.fastq \
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
--genomefile chr_A_B_unassigned.fasta --outdir ${outpath}Rust_smRNA_Predictions_ShortStack
```

##### Extract predicted siRNAs and miRNAs from ShortStack output files, get their composition and save their read counts for differential expression analysis later
```
cd Scripts
python ShortStack_to_Fasta.py ${outpath}Rust_smRNA_Predictions_ShortStack/Results.txt
python ShortStack_to_Fasta.py ${outpath}Rust_smRNA_Predictions_ShortStack/Results.txt

python Uracil_5Prime.py ${outpath}Rust_smRNA_Predictions_ShortStack/Results_siRNAs.fasta
python Uracil_5Prime.py ${outpath}Rust_smRNA_Predictions_ShortStack/Results_miRNAs.fasta

python Uracil_5Prime.py ${outpath}Wheat_smRNA_Predictions_ShortStack/Results_siRNAs.fasta
python Uracil_5Prime.py ${outpath}Wheat_smRNA_Predictions_ShortStack/Results_miRNAs.fasta

python ShortStack_to_Fasta_FilterCountMatrix.py ${outpath}Rust_smRNA_Predictions_ShortStack/Results.txt
python ShortStack_to_Fasta_FilterCountMatrix.py ${outpath}Wheat_smRNA_Predictions_ShortStack/Results.txt
cd ../
```

##### Differential expression analysis of rust and wheat sRNAs
```
EdgeR_Rust/Expression_smRNAs.R
EdgeR_Wheat/Expression_smRNAs.R
cd Scripts
python ShortStack_EdgeR_Rust.py ${outpath}
python ShortStack_EdgeR_Wheat.py ${outpath}
cd ../
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

##### Analysis of homologous sRNA loci in rust
```
cd Scripts
python Homologous_sRNAs.py smRNAs_upSpores.fasta
python Homologous_sRNAs.py smRNAs_up_early_infection.fasta
python Homologous_sRNAs.py smRNAs_up_late_infection.fasta
python Homologous_sRNAs.py smRNAs_noDE.fasta
cd ../
```

##### Get the chromosome locations for the sRNAs, so that they can be plotted onto the chromosomes (used in Figures 6 and 8)
```
cd Scripts
python KaryplotR_Repeats.py
python KaryplotR_Genes.py
python KaryplotR_smRNAs.py ${outpath}Rust_smRNA_PredictionsShortStack/smRNAs_up_early_infection.fasta
python KaryplotR_smRNAs.py ${outpath}Rust_smRNA_PredictionsShortStack/smRNAs_up_late_infection.fasta
python KaryplotR_smRNAs.py ${outpath}Rust_smRNA_PredictionsShortStack/smRNAs_upSpores.fasta
python KaryplotR_smRNAs.py ${outpath}Rust_smRNA_PredictionsShortStack/smRNAs_noDE.fasta
cd ../
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
