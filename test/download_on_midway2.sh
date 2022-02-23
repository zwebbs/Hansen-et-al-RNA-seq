# File Name: download_on_midway2.sh
# Created On: 2022-02-22
# Created By: ZW
# Purpose: runs the download script for test data for RNAseq pipeline
# should download six samples (3-treatment, 3-control) of human PE RNA-seq 

module load python/cpython-3.7.0
cd /project2/nobrega/zach/RNA_Seq_Pipeline/

python extra/SRA_fastq_download.py extra/acc_test.txt --prefetch-output-dir test/data/raw/ --fasterq-dump-output-dir test/data/raw/
