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
repeats="REPET/chrs_repet_no_SSR.bed"
genes="GeneAnnotation/Puccinia_graminis_tritici_21-0.sorted.bed"

# Set the overlap thresholds for smRNAs mapping to gene/repeats
a_overlap=0.33
# sRNAs overlapping/containing repeats/genes
#-f 	Minimum overlap required as a fraction of A
#-F 	Minimum overlap required as a fraction of B
declare -a StringArray=("smRNAs_upSpores" "smRNAs_up_early_infection" "smRNAs_up_late_infection" "smRNAs_noDE")
echo 
echo "---------------------"
for val in "${StringArray[@]}"; do
	echo 'smRNAs overlapping with repeats/genes'
	echo $val

	bedtools intersect -nonamecheck -f $a_overlap -a smRNA_Locations/${val}.sorted.bed -b ${repeats} -wo > smRNA_Locations/${val}_repeats_overlapping.bed
	bedtools intersect -nonamecheck -f $a_overlap -a smRNA_Locations/${val}.sorted.bed -b ${genes} -wo > smRNA_Locations/${val}_genes_overlapping.bed

done
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

##### 


