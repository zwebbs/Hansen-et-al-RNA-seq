# RNA_Seq_Pipeline

# Testing the pipeline on Gardner HPC Cluster

```
    -w <working directory for the analysis> [required]
    -c <analysis configuration file .json> [required]
    -d (BOOLEAN flag to complete a snakemake dry run) [optional]
```

## Running the Pipeline on Gardner HPC Cluster. 

Because we are unable to access the `qsub` utility on the compute nodes, I recommend running the snakemake manager process from your project directory on the head (login) node. This should require almost none of the head node's resources to manage, and should not run afoul of the CRI staff's preferences. Currently this is the only way I can think of to properly utilize the snakemake workflow management process on Gardner. 

```bash
./Pipeline_Submit_gardnerCRI.sh \
        -w `pwd`\
        -d \
        -c <config file name>.json \
        1> <stdout log file prefix>.out \
        2> <stderr log file prefix>.err
```
Conducting a dry run using the snakemake `--dry-run` feature allows you to check any syntax errors in the snakemake file or configuration file and is indicated by using the `-d` flag in the passed bash arguments.

```bash
./Pipeline_Submit_gardnerCRI.sh \
        -w `pwd`\
        -d \
        -c <config file name>.json \
        1> <stdout log file prefix>.out \
        2> <stderr log file prefix>.err
```
