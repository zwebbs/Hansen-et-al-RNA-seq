# RNA_Seq_Pipeline

# Testing the pipeline on Gardner HPC Cluster
Conducting a dry run using the snakemake `--dry-run` feature.
````
qsub Pipeline_Submit_gardnerCRI.sh -F "-w `pwd` -c config/RNA_Seq_Test.config.json -d"
````
