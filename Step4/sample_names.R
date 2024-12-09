library("dplyr")

# setwd
setwd("/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Code/exon_analysis/Data/feature_counts_output_exons")

# read in counts and sample names
sample_names <- read.csv("/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Code/exon_analysis/Data/AM_samples.csv")
count_data <- read.csv("mRNA_raw_reads.csv")

count_cols <- colnames(count_data)
count_cols <- count_cols[2:length(count_cols)]

new_names <- c()
for (col in count_cols) {
    new_name <- sample_names$sample_name[sample_names$alt_name==col]
    new_names <- c(new_names, new_name)
}

# rename count data cols
new_names <- c("Geneid", new_names)
colnames(count_data) <- new_names

# # reorder:
# ordered_cols <- sample_names$named_sample
# ordered_cols <- c("Geneid", ordered_cols)

# count_data <- count_data[, ordered_cols]

write.csv(count_data, "mRNA_raw_reads_named.csv", row.names=FALSE)
