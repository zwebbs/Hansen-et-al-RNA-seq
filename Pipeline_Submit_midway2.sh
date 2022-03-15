#!/bin/bash
#SBATCH --job-name=snakemake_RNA_Seq_Pipeline
#SBATCH --output=snakemake_RNA_Seq_Pipeline.out
#SBATCH --error=snakemake_RNA_Seq_Pipeline.err
#SBATCH --time=4:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G

# File Name:Pipeline_Submit_midway2.sh
# Created On: 2022-02-22
# Created By: ZW
# Purpose: runs the RNA-seq pipeline for given samples on midway2


# load required modules 
module load python/cpython-3.7.0
module load STAR/2.7.7a

# run the snakemake workflow
snakemake --snakefile RNA_Seq_Pipeline.snakefile \
    -j 28 -kp --rerun-incomplete --configfile test/RNA_Seq_PE_config_test.json \
    --cluster "sbatch \
     --nodes={resources.nodes} \
     --ntasks-per-node={resources.ntasks_per_node} \
     --cpus-per-task={resources.cpus_per_task} \
     --mem-per-cpu={resources.mem_per_cpu} \
     --job-name={rulename} \
     --error={rulename}.err \
     --output={rulename}.out \
     --time={resources.time} \
     --partition=broadwl" \
     $*



