# Human‑Read Extraction and Genetic Sex Inference from Gut Metagenomes

Bioinformatics workflow employed for the extraction of human-origin reads from shotgun metagenomic sequencing data

## Overview
![image](https://github.com/user-attachments/assets/9449c300-4b7f-4f61-b9e1-061d83d2b019)


This repository provides a **reproducible, end‑to‑end bioinformatics workflow** to

1. **Detect and remove human‑derived reads** from high‑throughput gut shotgun metagenomes;
2. **Quantify chromosome‑specific read depth**; and
3. **Infer genetic sex** using the **Y‑to‑X read‑depth ratio `(RYX)`**.

Although host DNA is routinely discarded as laboratory or computational
*contamination*, these reads can expose highly sensitive information (e.g.
biological sex, ancestry, pathogenic variants).  Our pipeline demonstrates that
`RYX` is a robust biomarker for sex prediction across different DNA‑extraction
chemistries and sequencing strategies (paired‑end and single‑end) while
highlighting the **ethical obligations** of handling human genetic material in
metagenomic studies.

---

## Biological & Ethical Background

* **Unintentional human DNA capture** is common in faecal metagenomic
  sequencing because epithelial cells and leukocytes are shed into stool.
* These reads can be exploited to infer individual traits (sex, ancestry,
  disease alleles) and therefore constitute **personally identifiable genomic
  information (PIGI)**.
* Bioinformatic de‑host pipelines are thus required **both** for improving
  microbial taxonomic accuracy **and** for safeguarding participant privacy.
* We standardise four widely used stool DNA‑extraction protocols and benchmark
  them under both PE and SE sequencing to assess residual host signal.

---

## Methodological Summary

**Step 1**: Build bowtie2 index
**Step 2**: Align metagenomic reads to Human reference genome
**Step 3**: Convert SAM to BAM and sorted the reads
**Step 4**: Extract aligned reads to human reference genome
**Step 5**: Convert BAM to FASTQ format and Compress the file as .gz
**Step 6**: Aligned compress read files withe the clad specific marker gene to check microbial contamination
**Step 7**: Marked duplicate reads with the .bam files 
**Step 8**: Extraction of chromosomal read count of the human readds


| Step | Tool(s)                           | Purpose                                            |
| ---- | --------------------------------- | -------------------------------------------------- |
| 1    | `fastqc`                          | Quality control of raw reads                       |
| 2    | `multiqc`                         | Aggregate QC reports                               |
| 3    | `bowtie2‑build`                   | Index GRCh38 or CHM13 reference genome             |
| 4    | `trimmomatic`, `fastp` (optional) | Adapter and quality trimming                       |
| 5    | `bowtie2`                         | Align reads to the human reference                 |
| 6    | `samtools view -F4`               | Isolate mapped (human) reads                       |
| 7    | `bedtools bamtofastq` + `pigz`    | Convert BAM → (paired) FASTQ for downstream checks |
| 8    | `metaphlan`                       | Check microbial contamination
| 9    | `picard MarkDuplicates`           | Mark PCR/optical duplicates                        |
| 10    | `samtools idxstats`               | Obtain per‑chromosome depth                        |
| 11   | `python (RYX) = readsY / readsX`  | Logistic model predicts sex                        |



---
![image](https://github.com/user-attachments/assets/905e7835-d1ec-4862-a90c-ea29f1a9b01f)


## Repository Structure

```
├── reference_genome/     # Index file of `human_reference_genome`
├── script/               # python script of sex prediction
├── data/                 # Data use to trained and test the model 
├── raw script/           # Raw script for extraction of human reads
├── workflow.sh           # One‑shot bash workflow (PE & SE)
└── README.md             # This file
```
## Usage Instructions for the Workflow Summary with command line
1. Build Bowtie2 index from GRCh38/chm13v2.0.fa genome
```shell
bowtie2-build GRCh38.primary_assembly.genome.fa human_genome_index_hg38
```  
2. Trim PE reads
```shell
for i in {19..21}; do [[ -f ${i}_R1.fastq.gz && -f ${i}_R2.fastq.gz ]] && trimmomatic PE -phred33 ${i}_R1.fastq.gz ${i}_R2.fastq.gz ${i}_R1_paired.fastq.gz ${i}_R1_unpaired.fastq.gz ${i}_R2_paired.fastq.gz ${i}_R2_unpaired.fastq.gz ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60; done
``` 
3. Trim SE reads
```shell
for i in {19..21}; do [[ -f ${i}.fastq.gz ]] && trimmomatic SE -phred33 ${i}.fastq.gz ${i}_trimmed.fastq.gz ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60; done
``` 
4. Align PE reads to human genome
```shell
for i in {02..54}; do [[ -f ${i}_R1.fastq.gz && -f ${i}_R2.fastq.gz ]] && bowtie2 --no-discordant -x chm13v2.0.fa/chm13v2.0_index -1 ${i}_R1.fastq.gz -2 ${i}_R2.fastq.gz -S ${i}_human.sam -p 16; done
```
5. Align SE reads to human genome
```shell 
for i in {02..54}; do [[ -f ${i}.fastq.gz ]] && bowtie2 -x chm13v2.0.fa/chm13v2.0_index -U ${i}.fastq.gz -S ${i}_human.sam -p 16; done 
```
6. Convert SAM to BAM, Sort BAM, and Index BAM 
```shell
for i in {14..32}; do [[ -f ${i}_human.sam ]] && samtools view -@ 16 -bS ${i}_human.sam > ${i}_human.bam; done
```
```shell 
for i in {14..32}; do [[ -f ${i}_human.bam ]] && samtools sort -@ 16 ${i}_human.bam -o ${i}_human_sorted.bam; done 
```
```shell 
for i in {14..32}; do [[ -f ${i}_human_sorted.bam ]] && samtools index -@ 16 ${i}_human_sorted.bam; done
```
7. Extract mapped reads only
```shell
for i in {01..16}; do [[ -f ${i}_human_sorted.bam ]] && samtools view -@ 16 -b -F 4 ${i}_human_sorted.bam > ${i}_human_reads.bam; done 
``` 
8. Convert BAM to gzipped paired FASTQ
```shell
for i in {14..32}; do [[ -f ${i}_human_reads.bam ]] && bedtools bamtofastq -i ${i}_human_reads.bam -fq >(pigz -p 16 > ${i}_new_R1.fastq.gz) -fq2 >(pigz -p 16 > ${i}_new_R2.fastq.gz); done
```
9. Check microbial contamination of paired FASTQ
```shell
for i in {14..32}; do [[ -f ${i}_R1.fastq.gz ]] && [[ -f ${i}_R2.fastq.gz ]] && metaphlan ${i}_R1.fastq.gz,${i}_R2.fastq.gz --bowtie2out ${i}.bowtie2.bz2 --nproc 20 --input_type fastq -o profiled_${i}.txt -t rel_ab_w_read_stats; done
```
10. Convert BAM to gzipped SE FASTQ
```shell
 for i in {02..54}; do [[ -f ${i}_human.bam ]] && bedtools bamtofastq -i ${i}_human.bam -fq >(pigz -p 4 > ${i}_human.fastq.gz); done
```
11. Check microbial contamination of SE FASTQ
```shell
for i in {14..32}; do [[ -f ${i}.fastq.gz ]] && metaphlan ${i}.fastq.gz --bowtie2out ${i}.bowtie2.bz2 --nproc 20 --input_type fastq -o profiled_${i}.txt -t rel_ab_w_read_stats; done
``` 
12. Remove  duplicates
```shell
for i in {47..54}; do [[ -f ${i}_human_reads.bam ]] && picard MarkDuplicates I=${i}_human_reads.bam O=${i}_deduplicated_human_reads.bam M=${i}_metrics.txt VALIDATION_STRINGENCY=LENIENT; done
```
13. Filter duplicates (flag 0x400) 
```shell
for i in {14..32}; do [[ -f ${i}_deduplicated_human_reads.bam ]] && samtools view -h -b -F 0x400 ${i}_deduplicated_human_reads.bam > ${i}_final_human_reads.bam; done
```  
14. Index final BAM
```shell
for i in {14..32}; do [[ -f ${i}_final_human_reads.bam ]] && samtools index -@ 16 ${i}_final_human_reads.bam; done
```  
15. Get chromosomal read counts
```shell
for i in {14..32}; do [[ -f ${i}_final_human_reads.bam ]] && samtools idxstats ${i}_final_human_reads.bam > ${i}_chromosomal_read_count.txt; done 
```
16. Get per-chromosome depth from BAM over Y regions
```shell
bedtools genomecov -ibam no_duplicates.bam -g Y.bed -bg | awk '{chr[$1]+=$4} END {for (c in chr) print c, chr[c]}'
```  

17. To extract chromosomal read count from the merge_chromosome.txt
```shell
for file in *.txt; do
    echo -e "$file" > temp_"$file"  
    awk '{print $3}' "$file" >> temp_"$file" 
done
```
```shell
for file in *.txt; do
    echo -e "$file" > temp_"$file"  
    awk '{print $3}' "$file" >> temp_"$file" 
done
```
```shell
paste temp_*.txt > merged_3rd_column.txt
```
```shell
rm temp_*.txt
```
```shell
paste -d ',' temp_*.txt > merged_3rd_column.txt
```
18. For all the column merge
```shell
for file in *.txt; do
    echo -e "$file" | cat - "$file" > temp_"$file"
done
```
```shell
paste temp_*.txt > merged_chromosomal_read_count.txt
```
```shell
rm temp_*.txt
```
```shell
paste -d ',' temp_*.txt > merged_chromosomal_read_count.txt
```

---

## Installation

```bash
conda env create -f env/host_removal.yml
conda activate host_removal
```

The environment installs **bowtie2 ≥2.5**, **samtools ≥1.20**, **bedtools ≥2.31**,
**picard ≥3.0**, **trimmomatic ≥0.40**, **pigz**, **python ≥3.11** with
`pandas`, `numpy`, `matplotlib`, `scikit‑learn`.

---
## Notes
In conclusion, this work demonstrates that the `Y : X` chromosomal read-depth ratio is a robust proxy for inferring genetic sex from gut metagenomic data. At the same time, it underscores the expanded ethical obligations inherent to high-resolution sequencing: the stewardship of human-derived reads must rigorously follow ethical and regulatory standards, both to uphold scientific integrity and to maintain public confidence in biomedical research.
## Contact Information
For questions or issues regarding this workflow, please contact Sahid Afrid Mollick (SAM) (sahidafridm@gmail.com).







