#!/bin/bash
#PBS -N snakemake_RNA_Seq_Pipeline
#PBS -S /bin/bash
#PBS -l walltime=08:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -o snakemake_RNA_Seq_Pipeline.out
#PBS -e snakemake_RNA_Seq_Pipeline.err

# File Name: Pipeline_Submit_gardnerCRI.sh
# Created On: 2022-03-29
# Created By: ZW
# Purpose: runs the RNA-seq pipeline for given samples on CRI Gardner HPC

# check passed commandline arguments
## of which there are three.
###   -w <working directory for the analysis> [required]
###   -c <analysis configuration file .json> [required]
###   -d (BOOLEAN flag to complete a snakemake dry run) [optional]

while getopts ":w:c:d" 'opt';
do
    case ${opt} in
        w)
            work_dir="${OPTARG}"
            ;; # capture the desired working directory and change to it.
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

# print the current working directory
cd ${work_dir} # change to desired work directory
echo "selected work dir: ${work_dir}"
echo "selected config file: ${config_file}"
echo "dry run ?: ${dry_run_flag}"

# load required modules
module load gcc/6.2.0
module load python/3.7.6

# run the snakemake workflow
snakemake --snakefile RNA_Seq_Pipeline.snakefile \
    -j 24 -kp --rerun-incomplete --configfile ${config_file} \
    --cluster "qsub -V -l walltime={resources.walltime} \
     -l nodes={resources.nodes}:ppn={resources.processors_per_node} \
     -l mem={resources.total_memory}mb \
     -N {rulename}_{resources.job_id} \
     -S /bin/bash \
     -e {resources.log_dir}{rulename}_{resources.job_id}.err \
     -o {resources.log_dir}{rulename}_{resources.job_id}.out" \
     ${dry_run_flag}
