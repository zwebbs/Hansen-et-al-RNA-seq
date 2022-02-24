#!/bin/bash
#SBATCH --job-name=SRA_download_fastq
#SBATCH --output=SRA_download_fastq.out
#SBATCH --error=SRA_download_fastq.err
#SBATCH --time=01:30:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=2000


# File Name: download_on_midway2.sh
# Created On: 2022-02-22
# Created By: ZW
# Purpose: runs the download script for test data for RNAseq pipeline
# should download six samples (3-treatment, 3-control) of human PE RNA-seq 

module load python/cpython-3.7.0
cd /project2/nobrega/zach/RNA_Seq_Pipeline/

python test/SRA_fastq_download.py test/acc_test2.txt --prefetch-output-dir test/data/raw/ --fasterq-dump-output-dir test/data/raw/
