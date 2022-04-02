#!/bin/bash
#PBS -N snakemake_RNA_Seq_Pipeline
#PBS -S /bin/bash
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -o snakemake_RNA_Seq_Pipeline.out
#PBS -e snakemake_RNA_Seq_Pipeline.err

# File Name: Pipeline_Submit_gardnerCRI.sh
# Created On: 2022-03-29
# Created By: ZW
# Purpose: runs the RNA-seq pipeline for given samples on CRI Gardner HPC

# check passed commandline arguments
## of which there are two. 
###   -c <analysis configuration file .json> [required]
###   -d (BOOLEAN flag to complete a snakemake dry run) [optional]

while getopts ":c:d" 'opt';
do
    case ${opt} in
        c) 
            config_file="${OPTARG}"
            ;; # capture command line argument for config file
        d) 
            dry_run_flag="--dry-run"
            ;; # capture boolean for whether or not its a 'dry-run'
        *) 
            echo "INVALID_OPTION -- ${OPTARG}"
            exit 1
            ;; # capture all other (invalid) inputs
    esac
done

# load required modules
module load gcc/6.2.0
module load python/3.7.6
module load STAR/2.6.1d
module load intel/2017
module load samtools/1.10
module load htslib/1.10.2

# run the snakemake workflow
snakemake --snakefile RNA_Seq_Pipeline.snakefile \
    -j 24 -kp --rerun-incomplete --configfile ${config_file} \
    --cluster "squb -l walltime={resources.walltime} \
     -l nodes={resources.nodes}:ppn={resources.processors_per_node} \
     -l mem={resources.total_memory}mb \
     -N {rulename}_{resources.job_id} \
     -S /bin/bash \
     -e {resources.logdir}{rulename}_{resources.job_id}.err \
     -o {resources.logdir}{rulename}_{resources.job_id}.out" \
     ${dry_run_flag}
