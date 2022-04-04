# RNA_Seq_Pipeline

# Testing the pipeline on Gardner HPC Cluster

```
    -w <working directory for the analysis> [required]
    -c <analysis configuration file .json> [required]
    -d (BOOLEAN flag to complete a snakemake dry run) [optional]
```

Conducting a dry run using the snakemake `--dry-run` feature using the `-d` flag in the passed bash arguments.
````bash
qsub Pipeline_Submit_gardnerCRI.sh \
         -F "-w `pwd` -c config/RNA_Seq_Test.config.json -d"
````
