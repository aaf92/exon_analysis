#!/bin/bash
#SBATCH --job-name=20241119_featureCounts_exons
#SBATCH --output=/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Data/AM/mRNA/feature_counts_output_exons/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --cluster=htc
#SBATCH --time=12:00:00

# load relevant modules
module purge
# module load gcc/10.2.0
# module load star/2.7.9a
# module load r/4.2.0
module load subread/2.0.1

#### Input:
gene_anno="/ix/djishnu/Aaron_F/TLR4_signaling/Code/scripts/Data_processing/STAR_aligner/mm10_build/mm10.refGene.gtf"
output_dir="/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Data/AM/mRNA/feature_counts_output_exons"
output_f="mRNA_raw_reads.txt"
fastq_dir="/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Data/AM/fastq"

# directory to feature counts output
mkdir -p "$output_dir"

featureCounts -T 8 \
-O \
-t exon \
-p \
-a $gene_anno \
-g exon_id \
-o "$output_dir/$output_f" \
"$fastq_dir"/*.bam