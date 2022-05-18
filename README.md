# RNA-Seq Pipeline

## Preparing the Workflow Environment

This Snakemake workflow requires Python 3.7, and was completed on Python 3.7.6. Later distributions of Python >= 3.7.6 may also work but we make no garuntees of forward compatibility. For maintaining multiple versions of Python on a single system, we highly recommend [Pyenv](https://github.com/pyenv/pyenv).  
With Python 3.7.6 as your primary interpreter, you can run the following commands in the terminal in order to create an isolated virtual environment and download the required python dependencies from the Python Package Index (PyPI).

```bash
-bash-4.1$ python3 -m venv env/
-bash-4.1$ source env/bin/activate
(env) -bash-4.1$ pip install -r requirements.txt
```

Additionally, we run these analyses with the help of a PBS-Torque scheduling system. Users not on a system where `qsub` style job submission is standard may have to adapt the workflow submission script and configuration resources to their own architecture.


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
