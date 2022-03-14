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
    -j 5 --configfile test/RNA_Seq_PE_config_test.json \
    --cluster-config test/RNA_Seq_PE_cluster_config_test.json \
    --cluster "sbatch \
     --nodes={cluster.nodes} \
     --ntasks-per-node={cluster.ntasks_per_node} \
     --cpus-per-task={cluster.cpu_per_tasks} \
     --mem-per-cpu={cluster.mem_per_cpu} \
     --job-name={cluster.job_name} \
     --error={cluster.error} \
     --output={cluster.output} \
     --time={cluster.time} \
     --partition=broadwl
    "



