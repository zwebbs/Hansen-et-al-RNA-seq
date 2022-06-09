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
## of which there are two.
###   -c <analysis configuration file .json> [required]
###   -d (BOOLEAN flag to complete a snakemake dry run) [optional]

while getopts ":j:c:d" 'opt';
do
    case ${opt} in
        j)
            parallel_jobs=${OPTARG}
            ;; # capture command line argument for number of concurrent jobs
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

# print options
echo "number of concurrent jobs: ${parallel_jobs}"
echo "selected config file: ${config_file}"
echo "dry run ?: ${dry_run_flag}"

# run the snakemake workflow
snakemake --snakefile Snakefile \
    -j ${parallel_jobs} -kp --rerun-incomplete \
    --config cluster_id="GARDNER" extended_config=${config_file} \
    --cluster "qsub -V -l walltime={resources.walltime} \
     -l nodes={resources.nodes}:ppn={resources.processors_per_node} \
     -l mem={resources.total_memory}mb \
     -N {rulename}_{resources.job_id} \
     -S /bin/bash \
     -e {resources.log_dir}{rulename}_{resources.job_id}.err \
     -o {resources.log_dir}{rulename}_{resources.job_id}.out" \
     ${dry_run_flag} 1> "hansen-rnaseq_gardner.out" 2> "hansen-rnaseq_gardner.err"

