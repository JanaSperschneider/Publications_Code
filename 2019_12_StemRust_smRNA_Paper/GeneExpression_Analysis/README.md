#### RNAseq analysis of PGT 21-0 gene expression data

######  Run cutadapt and trimmomatic on the RNAseq reads

```
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW0U1.fastq.gz Rust_RNAseq/RW0U1_C8006ANXX_ATCACGAT_L001_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW0U2.fastq.gz Rust_RNAseq/RW0U2_C8006ANXX_CGATGTAT_L001_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW0U3.fastq.gz Rust_RNAseq/RW0U3_C8006ANXX_TTAGGCAT_L001_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW2I1.fastq.gz Rust_RNAseq/RW2I1_C8006ANXX_TGACCA_L002_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW2I2.fastq.gz Rust_RNAseq/RW2I2_C8006ANXX_ACAGTG_L002_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW2I3.fastq.gz Rust_RNAseq/RW2I3_C8006ANXX_GCCAAT_L003_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW3I1.fastq.gz Rust_RNAseq/RW3I1_C8006ANXX_CAGATC_L003_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW3I2.fastq.gz Rust_RNAseq/RW3I2_C8006ANXX_ACTTGA_L004_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW3I3.fastq.gz Rust_RNAseq/RW3I3_C8006ANXX_GATCAG_L004_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW4I1.fastq.gz Rust_RNAseq/RW4I1_C8006ANXX_TAGCTT_L005_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW4I2.fastq.gz Rust_RNAseq/RW4I2_C8006ANXX_GGCTAC_L005_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW4I3.fastq.gz Rust_RNAseq/RW4I3_C8006ANXX_CTTGTA_L005_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW5I1.fastq.gz Rust_RNAseq/RW5I1_C8006ANXX_AGTCAACA_L006_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW5I2.fastq.gz Rust_RNAseq/RW5I2_C8006ANXX_AGTTCCGT_L006_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW5I3.fastq.gz Rust_RNAseq/RW5I3_C8006ANXX_ATGTCAGA_L006_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW6I1.fastq.gz Rust_RNAseq/RW6I1_C8006ANXX_CCGTCCCG_L001_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW6I2.fastq.gz Rust_RNAseq/RW6I2_C8006ANXX_GTCCGCAC_L001_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW6I3.fastq.gz Rust_RNAseq/RW6I3_C8006ANXX_GTGAAACG_L001_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW7I1.fastq.gz Rust_RNAseq/RW7I1_C8006ANXX_GTGGCCTT_L001_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW7I2.fastq.gz Rust_RNAseq/RW7I2_C8006ANXX_GTTTCGGA_L001_R1.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 25 -o Rust_RNAseq/cutadapt/RW7I3.fastq.gz Rust_RNAseq/RW7I3_C8006ANXX_CGTACGTA_L001_R1.fastq.gz

trim_galore --paired Rust_RNAseq/Sample_Nmu_210_HSRNA_1/Nmu_210_HSRNA_1_CCGTCC_L005_R1_001.fastq.gz Rust_RNAseq/Sample_Nmu_210_HSRNA_1/Nmu_210_HSRNA_1_CCGTCC_L005_R2_001.fastq.gz -o Rust_RNAseq/trimgalore
trim_galore --paired Rust_RNAseq/Sample_Nmu_210_HSRNA_2/Nmu_210_HSRNA_2_GTCCGC_L005_R1_001.fastq.gz Rust_RNAseq/Sample_Nmu_210_HSRNA_2/Nmu_210_HSRNA_2_GTCCGC_L005_R2_001.fastq.gz -o Rust_RNAseq/trimgalore
trim_galore --paired Rust_RNAseq/Sample_Nmu_210_HSRNA_3/Nmu_210_HSRNA_3_GTGAAA_L005_R1_001.fastq.gz Rust_RNAseq/Sample_Nmu_210_HSRNA_3/Nmu_210_HSRNA_3_GTGAAA_L005_R2_001.fastq.gz -o Rust_RNAseq/trimgalore

trim_galore --paired Rust_RNAseq/Sample_Nmu_PGTGS_1/NMu_PGTGS_1_AGTCAA_L005_R1_001.fastq.gz Rust_RNAseq/Sample_Nmu_PGTGS_1/NMu_PGTGS_1_AGTCAA_L005_R2_001.fastq.gz -o Rust_RNAseq/trimgalore
trim_galore --paired Rust_RNAseq/Sample_Nmu_PGTGS_2/Nmu_PGTGS_2_AGTTCC_L005_R1_001.fastq.gz Rust_RNAseq/Sample_Nmu_PGTGS_2/Nmu_PGTGS_2_AGTTCC_L005_R2_001.fastq.gz -o Rust_RNAseq/trimgalore
trim_galore --paired Rust_RNAseq/Sample_Nmu_PGTGS_3/Nmu_PGTGS_3_ATGTCA_L005_R1_001.fastq.gz Rust_RNAseq/Sample_Nmu_PGTGS_3/Nmu_PGTGS_3_ATGTCA_L005_R2_001.fastq.gz -o Rust_RNAseq/trimgalore
```

###### Align against transcripts with Salmon

```
salmon index -t Puccinia_graminis_tritici_21-0.transcripts.fa -i SalmonIndex/PGT_210_AGO_Fixed

for fn in Rust_RNAseq/cutadapt/*.fastq.gz
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant --validateMappings -i SalmonIndex/PGT_210_AGO_Fixed --libType SR -r ${fn} -p 4 -o quants/${samp}_quant --writeMappings=quants/${samp}_quant/mapping.sam 
done

f1=Rust_RNAseq/trimgalore/Nmu_210_HSRNA_1_CCGTCC_L005_R1_001_val_1.fq.gz
f2=Rust_RNAseq/trimgalore/Nmu_210_HSRNA_1_CCGTCC_L005_R2_001_val_2.fq.gz
salmon quant --validateMappings -i SalmonIndex/PGT_210_AGO_Fixed -l A -1 $f1 -2 $f2 -p 4 -o quants/HSRNA1_quant --writeMappings=quants/HSRNA1_quant/mapping.sam

f1=Rust_RNAseq/trimgalore/Nmu_210_HSRNA_2_GTCCGC_L005_R1_001_val_1.fq.gz
f2=Rust_RNAseq/trimgalore/Nmu_210_HSRNA_2_GTCCGC_L005_R2_001_val_2.fq.gz
salmon quant --validateMappings -i SalmonIndex/PGT_210_AGO_Fixed -l A -1 $f1 -2 $f2 -p 4 -o quants/HSRNA2_quant --writeMappings=quants/HSRNA2_quant/mapping.sam 

f1=Rust_RNAseq/trimgalore/Nmu_210_HSRNA_3_GTGAAA_L005_R1_001_val_1.fq.gz
f2=Rust_RNAseq/trimgalore/Nmu_210_HSRNA_3_GTGAAA_L005_R2_001_val_2.fq.gz
salmon quant --validateMappings -i SalmonIndex/PGT_210_AGO_Fixed -l A -1 $f1 -2 $f2 -p 4 -o quants/HSRNA3_quant --writeMappings=quants/HSRNA3_quant/mapping.sam

f1=Rust_RNAseq/trimgalore/NMu_PGTGS_1_AGTCAA_L005_R1_001_val_1.fq.gz
f2=Rust_RNAseq/trimgalore/NMu_PGTGS_1_AGTCAA_L005_R2_001_val_2.fq.gz
salmon quant --validateMappings -i SalmonIndex/PGT_210_AGO_Fixed -l A -1 $f1 -2 $f2 -p 4 -o quants/PGTGS1_quant --writeMappings=quants/PGTGS1_quant/mapping.sam

f1=Rust_RNAseq/trimgalore/Nmu_PGTGS_2_AGTTCC_L005_R1_001_val_1.fq.gz
f2=Eust_RNAseq/trimgalore/Nmu_PGTGS_2_AGTTCC_L005_R2_001_val_2.fq.gz
salmon quant --validateMappings -i SalmonIndex/PGT_210_AGO_Fixed -l A -1 $f1 -2 $f2 -p 4 -o quants/PGTGS2_quant --writeMappings=quants/PGTGS2_quant/mapping.sam

f1=Rust_RNAseq/trimgalore/Nmu_PGTGS_3_ATGTCA_L005_R1_001_val_1.fq.gz
f2=Rust_RNAseq/trimgalore/Nmu_PGTGS_3_ATGTCA_L005_R2_001_val_2.fq.gz
salmon quant --validateMappings -i SalmonIndex/PGT_210_AGO_Fixed -l A -1 $f1 -2 $f2 -p 4 -o quants/PGTGS3_quant --writeMappings=quants/PGTGS3_quant/mapping.sam

```

Salmon_DE_Analysis.R
