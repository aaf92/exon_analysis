# Exon Analysis Workflow

This document outlines the steps for analyzing exon-level RNA-seq data, including alignment, feature counting, normalization, and visualization.

---

## **1. Alignment**
- **Script:** `star_align.slurm`  
- **Input:** FASTQ files  
- **Output:** BAM files  

---

## **2. Exon-Level Feature Counts**
- **Script:** `feature_counts_exons.slurm`  
- **Input:** `*.bam`  
- **Output:** `mRNA_raw_reads.txt`  

---

## **3. Convert Feature Counts to CSV**
- **Script:** `get_subread_countTable.rmd`  
- **Input:** `mRNA_raw_reads.txt`  
- **Output:** `mRNA_raw_reads.csv`  

---

## **4. Map FASTQ File Names to Sample Names**
- **Script:** `sample_names.R`  
- **Input:**  
  - `mRNA_raw_reads.csv`  
  - `AM_samples.csv`  
- **Output:** `mRNA_raw_reads_named.csv`  

---

## **5. Split Samples into Respective Experiments**
- **Script:** `separate_experiments.R`  
- **Input:** `mRNA_raw_reads_named.csv`  
- **Output:**  
  - `mRNA_raw_reads_named_exp1.csv`  
  - `mRNA_raw_reads_named_exp2.csv`
- **Notes:**  
  - Modify exp_num variable 

---

## **6. Remove Low-Quality Samples**
- **Script:** `remove_samples.R`  
- **Input:**  
  - `mRNA_raw_reads_named_exp{1 or 2}.csv`  
  - **Samples to Remove:**  
    - **Experiment 1:** `IL33_WT_24_Exp1_Rep3`  
    - **Experiment 2:** `Staph_GATA2_24_Exp2_Rep2`  
- **Output:**  
  - `mRNA_raw_reads_named_exp1_removedSamples.csv`  
  - `mRNA_raw_reads_named_exp2_removedSamples.csv`
- **Notes:**  
  - Modify exp_num variable 

---

## **7. Normalize Each Experiment Using DESeq2**
- **Script:** `normalize_samples.rmd`  
- **Input:**  
  - `mRNA_raw_reads_named_exp{1 or 2}_removedSamples.csv`  
- **Output:**  
  - `{YYYYMMDD}_deseq2_normalized_RNAseq_counts_exp1_removedSamples.csv`  
  - `{YYYYMMDD}_deseq2_normalized_RNAseq_counts_exp2_removedSamples.csv`
- **Notes:**  
  - Modify date and exp_num variable 

---

## **8. Generate Heatmap for Exons of Interest**
- **Script:** `heatmap_exons.ipynb`  
- **Input:**  
  - `{YYYYMMDD}_deseq2_normalized_RNAseq_counts_exp{1 or 2}_removedSamples.csv`  
  - `exons_to_plot.txt`  
  - `gata2_exon_names.csv`  
- **Output:**  
  - `gata2_exon_heatmap_deseq2_normalized.pdf`  
- **Notes:**  
  - Make sure to modify the inputs cell variables: exp_num, results_dir, data_dir
