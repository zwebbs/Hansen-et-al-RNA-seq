#!/bin/bash
#SBATCH --job-name=RNA_Seq_Pipeline_Test
#SBATCH --output=RNA_Seq_Pipeline_Test.out
#SBATCH --error=RNA_Seq_Pipeline_Test.err
#SBATCH --time=02:30:00
#SBATCH --partition=broadwl
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000


# File Name:Pipeline_Submit_midway2.sh
# Created On: 2022-02-22
# Created By: ZW
# Purpose: runs the RNA-seq pipeline for given samples on midway2


# load required modules 
module load python/cpython-3.7.0
module load STAR/2.6.1b

# run the snakemake workflow
snakemake --snakefile RNA_Seq_Pipeline.snakefile \
    -j 5 --configfile test/RNA_Seq_PE_config_test.json



