TLR4_signaling/ATAC_data/SRR_download/samtools/bam_index.slurm#!/bin/bash
#SBATCH --job-name=20240722_BMDM_fastqAlign
#SBATCH --output=/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Data/BMDM/fastq/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --cluster=htc
#SBATCH --time=12:00:00

# load relevant modules
module purge
# module load gcc/10.2.0
module load star/2.7.9a
# module load r/4.2.0

#### Inputs:
# specify if paired as true or false:
paired=true

# specify if fastq files are zipped or not
gzipped=true

# directories containing fastq files and the genome for alignment
fastq_dir="/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Data/BMDM/fastq"
genome_dir="/ix/djishnu/Aaron_F/morgan_satarupa_RNAseq/Code/scripts/Data_processing/STAR_aligner/mm10_build/mm10_index"

#### End Inputs

#### Main Script
# go to fastq directory
cd $fastq_dir

# get array of fastq files
if [ "$gzipped" = true ]; then
    if [ "$paired" = true ]; then
        # Create arrays for read1 and read2 files
        read1_files=(*R1*.fastq.gz)
        read2_files=(*R2*.fastq.gz)
    else
        read1_files=(*_pass.fastq.gz)
    fi
else
    if [ "$paired" = true ]; then
        # Create arrays for read1 and read2 files
        read1_files=(*R1*.fastq)
        read2_files=(*R2*.fastq)
    else
        read1_files=(*_pass.fastq)
    fi
fi


# iterate through fastq and run star align
if [ "$paired" = true ]; then
    # Sort arrays naturally
    IFS=$'\n' read1_files_sorted=($(sort -V <<<"${read1_files[*]}"))
    IFS=$'\n' read2_files_sorted=($(sort -V <<<"${read2_files[*]}"))

    # Iterate over the arrays
    for ((i=0; i<${#read1_files_sorted[@]}; i++)); do
        read1_file="${read1_files_sorted[i]}"
        read2_file="${read2_files_sorted[i]}"

        if [ "$gzipped" = true ]; then
            # Run STAR alignment for this pair of read files
            STAR --genomeDir $genome_dir \
            --runThreadN 8 \
            --readFilesIn "$read1_file" \
            --readFilesCommand gunzip -c \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix $fastq_dir/${read1_file}_
        else
            # Run STAR alignment for this pair of read files
            STAR --genomeDir $genome_dir \
            --runThreadN 8 \
            --readFilesIn "$read1_file" \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix $fastq_dir/${read1_file}_
        fi
    done
else
    # Sort arrays naturally
    IFS=$'\n' read1_files_sorted=($(sort -V <<<"${read1_files[*]}"))

    # Iterate over the arrays
    for ((i=0; i<${#read1_files_sorted[@]}; i++)); do
        read1_file="${read1_files_sorted[i]}"

        if [ "$gzipped" = true ]; then
            # Run STAR alignment for this pair of read files
            STAR --genomeDir $genome_dir \
            --runThreadN 8 \
            --readFilesIn "$read1_file" \
            --readFilesCommand gunzip -c \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix $fastq_dir/${read1_file}_
        else
            # Run STAR alignment for this pair of read files
            STAR --genomeDir $genome_dir \
            --runThreadN 8 \
            --readFilesIn "$read1_file" \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix $fastq_dir/${read1_file}_
        fi
    done
fi
