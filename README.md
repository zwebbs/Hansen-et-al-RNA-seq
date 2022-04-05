# RNA-Seq Pipeline


## Running the Pipeline 

## Running the Pipeline on Gardner HPC Cluster

To run the Pipeline on the Gardner Cluster, we have supplied a job manager scipt `Pipeline_Submit_gardnerCRI.sh`. Executing this script from the project directory should be all that is required --assuming that the configuration has been properly completed. The commandline options for `Pipeline_Submit_gardnerCRI.sh` are listed below:
```
    -w <working directory for the analysis> [required]
    -c <analysis configuration file .json> [required]
    -d (BOOLEAN flag to complete a snakemake dry run) [optional]
```

Gardner makes use of the **qsub-Torque** job scheduler system, which is esily supported by Snakemake's `--cluster` management option. Because of specific restrictions on Gardner, we are unable to access the `qsub` utility on the compute nodes, thus I recommend running the snakemake manager script from your project directory on the head (login) node.
Using `qsub` to submit the job manager script will most likely result in a Error Code  (Exit 127) for a failed jobscript submission. Running it instead from the head node should require almost none of the head node's resources --thus not running afoul of the CRI staff's preferences. Currently this is the only way I can think of to properly utilize the snakemake workflow management process on Gardner. 

```bash
./Pipeline_Submit_gardnerCRI.sh \
        -w `pwd`\
        -c <config file name>.json \
        1> <stdout log file prefix>.out \
        2> <stderr log file prefix>.err
```
Conducting a dry run using the snakemake `--dry-run` feature allows you to check any syntax errors in the snakemake file or configuration file and is indicated by using the `-d` flag in the passed bash arguments.

```bash
./Pipeline_Submit_gardnerCRI.sh \
        -w `pwd`\
        -c <config file name>.json \
        -d \
        1> <stdout log file prefix>.out \
        2> <stderr log file prefix>.err
```
