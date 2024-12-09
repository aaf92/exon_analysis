library("dplyr")

# input:
dir <- "/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Code/exon_analysis/Data/feature_counts_output_exons"
exp_num <- "1"

# setwd
setwd(dir)

# read in counts
count_data <- read.csv("mRNA_raw_reads_named.csv")

# current cols
cur_cols <- colnames(count_data)

# subset to exp# cols
exp_cols <- grepl(paste0("Exp", exp_num), cur_cols)
exp_cols <- cur_cols[exp_cols]

new_cols <- c(cur_cols[1], exp_cols)

# subset data
subset_counts <- count_data[,new_cols]

write.csv(subset_counts, paste0("mRNA_raw_reads_named_exp", exp_num, ".csv"), row.names=FALSE)
