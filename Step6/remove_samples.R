library(dplyr)

#### inputs ####
# input:
dir <- "/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Code/exon_analysis/Data/feature_counts_output_exons"
exp_num <- "1"

#### setwd for output ease
setwd("/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Code/exon_analysis/Data/feature_counts_output_exons")

if (exp_num=="1") {
    counts_f <- "mRNA_raw_reads_named_exp1.csv"
    counts_out <- "mRNA_raw_reads_named_exp1_removedSamples.csv"
    #### samples to remove
    rem_samples <- c("IL33_WT_24_Exp1_Rep3")
}
else{
    counts_f <- "mRNA_raw_reads_named_exp2.csv"
    counts_out <- "mRNA_raw_reads_named_exp2_removedSamples.csv"
    #### samples to remove
    rem_samples <- c("Staph_GATA2_24_Exp2_Rep2")
}

#### main ####
#### read in counts
counts <- read.csv(counts_f)

#### remove samples from counts
counts_new <- counts %>% select(-all_of(rem_samples))

#### save
write.csv(counts_new, counts_out, row.names=FALSE)
