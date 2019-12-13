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
