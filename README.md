# Somatic Variants Pipeline

---

## Table of Contents

1. [Server Information](#1-server-information)
2. [Genomic Sequencing Data](#2-genomic-sequencing-data)
3. [Preprocessing Steps](#3-preprocessing-steps)
4. [Somatic Variant Calling and Filtering](#4-somatic-variant-calling-and-filtering)
5. [Annotate Variants](#5-annotate-variants)
6. [Copy Number Variation](#6-copy-number-variation)


---


## 1. Server Information

### Prerequisites
- **Tools:** Install  [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html), [FileZilla](https://filezilla-project.org/download.php), [IGV](https://igv.org/doc/desktop/)  
- **Server Details:**  
  - **IP Address:** `168.105.161.70`  
  - **Port:** `22`  
  - **Access:** Requires JABSOM or UH network (use VPN for remote access).
  - Note: With PuTTY and FileZilla you can connect to server.
 

### Security Practices
- Avoid multiple failed login attempts to prevent account locking.  
- Use strong passwords or SSH keys.  
- Log out after completing tasks.

### Installing Miniconda
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
### Installing software
```bash
# conda install mamba
# mamba install sra-tools fastqc trimmomatic multiqc curl 
# mamba install bwa samtools picard gatk4
# mamba install snpeff snpsift
```

---

## 2. Genomic Sequencing Data
Explore somatic mutations on HCC1143 tumor sample with matched normal sample.
It is the breast cancer cell line.

Starting files:
BAM files: 
/home/bqhs/mutect2/tumor.bam, 
/home/bqhs/mutect2/normal.bam

Sample names: HCC1143_tumor, HCC1143_normal

Reference: /home/bqhs/mutect2/Homo_sapiens_assembly38.fasta 
PoN: /home/bqhs/mutect2/chr17_m2pon.vcf.gz
Germline AF resource: /home/bqhs/mutect2/chr17_af-only-gnomad_grch38.vcf.gz
Intervals: /home/bqhs/mutect2/targets_chr17.interval_list
In the interest of time, we’ll only do somatic variant calling on parts of chr17
Common biallelic SNPs: /home/bqhs/mutect2/chr17_small_exac_common_3_grch38.vcf.gz (for GetPileupSummaries)

Run the pipeline: Mutect2, filter contamination and orientation bias
When running Mutect2, specify –L /home/bqhs/mutect2/chr17plus.interval_list
Examine results in IGV. Load:
The BAM files
Germline AF resource and PoN VCF
Your filtered Mutect2 results (VCF)
What do you observe at these sites?
chr17:7,666,402-7,689,550
chr17: 7,221,420; 19,748,387; 50,124,771

CNVKit
Starting files:
Use the same BAM files and reference genome as Mutect2 exercise
Target list: /home/bqhs/mutect2/targets_chr17.interval_list
Gene annotation: /home/bqhs/mutect2/refFlat.txt
Run CNVkit
Observe the effect of changing the purity setting when converting to copy number (integers)

---

## 3. Preprocessing Steps 

### Common to Germline and Somatic Variants
1. QC raw reads (FastQC, Trimmomatic) → Ensures high-quality reads
2. Align reads (bwa mem) → Maps reads to the reference genome
3. Sort alignment (samtools sort) → Organizes BAM file
4. Mark duplicates (picard MarkDuplicates) → Removes PCR artifacts
5. Recalibrate base quality scores (gatk BaseRecalibrator, gatk ApplyBQSR) → Adjusts systematic sequencing errors
6. Index BAM file (samtools index) → Enables efficient querying

---

## 4. Somatic Variant Calling and Filtering 

```
mkdir -p Somatic_mutation/Exercise
cd Somatic_mutation/Exercise
```

### 1. Run gatk Mutect2
+ Detects tumor-specific mutations
+ Requires both tumor and optionally normal BAM files for filtering

```
 gatk Mutect2 -R /home/bqhs/mutect2/Homo_sapiens_assembly38.fasta -I /home/bqhs/mutect2/tumor.bam -tumor HCC1143_tumor -I /home/bqhs/mutect2/normal.bam -normal HCC1143_normal -pon /home/bqhs/mutect2/chr17_m2pon.vcf.gz --germline-resource /home/bqhs/mutect2/chr17_af-only-gnomad_grch38.vcf.gz -L /home/bqhs/mutect2/chr17plus.interval_list -O somatic_m2.vcf.gz
```

```
# Sample name
samtools view -H /home/bqhs/mutect2/tumor.bam | grep '@RG'
samtools view -H /home/bqhs/mutect2/normal.bam | grep '@RG'
gatk GetSampleName -I /home/bqhs/mutect2/tumor.bam  -O tumor.txt
gatk GetSampleName -I /home/bqhs/mutect2/normal.bam  -O normal.txt
```

```
# subsets records that contain a comma in the 5th column.
zcat somatic_m2.vcf.gz | awk '$5 ~","'
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele-

zcat somatic_m2.vcf.gz | head -n 10000

# Can view the standard Format in the VCF header
zcat somatic_m2.vcf.gz | grep '##FORMAT'

# Can view the standard INFO in the VCF header
zcat somatic_m2.vcf.gz | grep '##INFO'
```

### 2. Run gatk GetPileupSummaries
+ Summarize pileup metrics 
+ (gatk GetPileupSummaries) → Gathers read count metrics

```
# For Normal
gatk GetPileupSummaries -R /home/bqhs/mutect2/Homo_sapiens_assembly38.fasta -I /home/bqhs/mutect2/normal.bam -V /home/bqhs/mutect2/chr17_small_exac_common_3_grch38.vcf.gz -L /home/bqhs/mutect2/targets_chr17.interval_list -O normal.pileups.table
```

```
cat tumor.pileups.table | grep '^chr17' | awk '$5>=3'
cat normal.pileups.table | grep '^chr17' | awk '$5>=3'
```


### 3. gatk CalculateContamination
+ Estimate sample contamination 
+ (gatk CalculateContamination) → Accounts for tumor impurity and normal contamination

```
 gatk CalculateContamination -I tumor.pileups.table -matched normal.pileups.table -O contamination.table
```

### 4. gatk FilterMutectCalls
+ Filter somatic variants 
(gatk FilterMutectCalls) → Removes low-confidence calls

```
gatk FilterMutectCalls -V somatic_m2.vcf.gz --contamination-table contamination.table -R /home/bqhs/mutect2/Homo_sapiens_assembly38.fasta -O somatic.filtered.vcf.gz
```

```
# contamination-table: Tables containing contamination information
zcat somatic.filtered.vcf.gz | grep '##FILTER'
zcat somatic.filtered.vcf.gz | grep 'germ-line'
zcat somatic.filtered.vcf.gz | grep 'PASS'
zcat somatic.filtered.vcf.gz | grep 'contamiation'
```

### 5. gatk CollectSequencingArtifactMetrics
+ Assess sequencing artifacts 
+ (gatk CollectSequencingArtifactMetrics) → Identifies systematic errors

```
gatk CollectSequencingArtifactMetrics -R /home/bqhs/mutect2/Homo_sapiens_assembly38.fasta -I /home/bqhs/mutect2/tumor.bam -EXT ".txt" -O tumor_artifact

```


---

## 5. Annotate Variants
### Functional annotation 
(snpEff) → Predicts functional effects of mutations

```
gunzip -c somatic.filtered.vcf.gz > somatic.filtered.vcf
snpEff -Xmx2g ann hg38 -v -s snpeff.html somatic.filtered.vcf > somatic.filtered.ann.vcf
SnpSift annotate /home/bqhs/hg38/dbsnp_146.hg38.vcf.gz somatic.filtered.ann.vcf > somatic.filtered_1.ann.vcf

bcftools query -l somatic.filtered_1.ann.vcf
```

### Extract relevant fields 
(SnpSift extractFields) → Filters and formats VCF for downstream analysis
```
SnpSift extractFields somatic.filtered_1.ann.vcf \
    ID CHROM POS REF ALT QUAL DP FILTER \
    ANN[0].GENE ANN[0].GENEID ANN[0].EFFECT ANN[0].IMPACT \
    ANN[0].BIOTYPE ANN[0].HGVS_C ANN[0].HGVS_P \
    GEN[0].GT GEN[0].GQ GEN[0].FT \
    GEN[1].GT GEN[1].GQ GEN[1].FT > somatic_final.txt
```
---

## 6. Copy Number Variation
### Run CNV detection 
(cnvkit.py batch) → Detects large deletions/amplifications

### Call CNVs 
(cnvkit.py call) → Determines CNV events per sample

---


---
