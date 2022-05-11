# File Name: Snakefile
# Created By: ZW
# Created On: 2022-02-17
# Purpose: Manages the workflow for RNA-seq alignment and quantification for
# both single-end and paired-end sequencing


# library + module imports
# ----------------------------------------------------------------------------
import gardnersnake as gs
from gardnersnake.misc.exceptions import UserError
from pathlib import Path

# workflow setup
# ----------------------------------------------------------------------------

# get appropriate configuration and metadata schemas according to which
# cluster we're using.
cluster_id = config["cluster_id"]
if cluster_id == "GARDNER":
    cfg_schema="CFG_GARDNER_BASIC"
    data_schema="META_GARDNER_SEQ_BASIC"
elif cluster_id == "MIDWAY2":
    cfg_schema = None  # TODO: implement gs midway2 schema and replace
    data_schema = None  # TODO: implement gs midway2 schema and replace
else:
    msg = f"Could not interpret cluster_id! - {cluster_id} not understood"
    raise UserError(msg)

# build configuration helper and metadata manager
cfg, meta = gs.read_yaml_extended_config(config['extended_config'])
CH = gs.ConfigurationHelper(cfg, schema_type=cfg_schema)
DM = gs.DataManager(meta, schema_type=data_schema)

# set analysis_id and working directory
analysis_id = CH.get_glob('analysis_id')
workdir: CH.get_glob('workdir')

# Globals
# --------------------------
ALIGN_DIR = CH.get_parameters("STAR_Align_Reads")["outdir_fprefix"]

print(DM.shared_data)
# workflow
#----------------

# Rule 0: Defining global outputs
# ----------------------------------------------------------------------------

rule all: # static output unpacking should all go last
    input:
        expand(ALIGN_DIR + "{run_id}.Aligned.sortedByCoord.out.bam", run_id=DM.get_shared_data({}, "run_name")),
        **DM.get_rule_data("STAR_Create_Genome_Index",["static_outputs"])


# Rule 1: Verify Genome Index or Create a new one.
# ----------------------------------------------------------------------------

# check if the index already exists and verify its contents
# use parameters and outputs of rule STAR_Create_Genome_Index because they
# are identical
genome_index_dir = CH.get_parameters("STAR_Create_Genome_Index")["genome_index_dir"]
if Path(genome_index_dir).exists():
    rule Verify_Index_Contents:
        params: **CH.get_parameters("STAR_Create_Genome_Index")
        output: **DM.get_rule_data("STAR_Create_Genome_Index", ["static_outputs"])
        resources: **CH.get_resources("Verify_Index_Contents")
        shell:
            "check_directory -o {output.rc_out}"
            " {params.genome_index_manifest} {params.genome_index_dir}"

else: # create the index if it doesn't exist, verify contents afterwords
    rule STAR_Create_Genome_Index:
        input: **DM.get_rule_data("STAR_Create_Genome_Index", ["static_inputs"])
        params: **CH.get_parameters("STAR_Create_Genome_Index")
        output: **DM.get_rule_data("STAR_Create_Genome_Index", ["static_outputs"])
        resources: **CH.get_resources("STAR_Create_Genome_Index")
        shell:
            "mkdir -p {params.genome_index_dir} && "
            "STAR --runThreadN {params.nthreads}"
            " --runMode genomeGenerate"
            " --genomeDir {params.genome_index_dir}"
            " --genomeFastaFiles {input.genome_fasta_path}"   
            " --sjdbGTFfile {input.transcript_gtf_path}"
            " --sjdbOverhang {params.sjdbOverhang}"
            " {params.extra_args} &&"
            " check_directory -o {output.rc_out}"
            " {params.genome_index_manifest}"
            " {params.genome_index_dir}"


# Rule 2: Aligning RNA-Seq reads using STAR
#-----------------------------------------------------------------------------

# scatter individual runs and align
rule STAR_Align_Reads:
    input:
        fastq1 = lambda wildcards: DM.get_shared_data({'run_name': f"{wildcards.run_id}"},"fastq1"),
        fastq2 = lambda wildcards: DM.get_shared_data({'run_name': f"{wildcards.run_id}"},"fastq2"),
        **DM.get_rule_data("STAR_Create_Genome_Index", ["static_outputs"])
    params:
        **CH.get_parameters("STAR_Align_Reads"),
        run_id = (lambda wildcards: wildcards.run_id),
        genome_index_dir=CH.get_parameters("STAR_Create_Genome_Index")["genome_index_dir"]
    resources:
        ** CH.get_resources("STAR_Align_Reads", return_job_id=False),
        job_id = (lambda wildcards: wildcards.run_id)
    output: CH.get_parameters("STAR_Align_Reads")["outdir_fprefix"] + "{run_id}.Aligned.sortedByCoord.out.bam"
    shell:
        "STAR --runThreadN {params.nthreads}"
        " --runMode {params.run_mode}"
        " --quantMode {params.quant_mode}"
        " --genomeDir {params.genome_index_dir}/"
        " --readFilesIn {input.fastq1} {input.fastq2}"
        " --readFilesCommand {params.read_files_command}"
        " --outFileNamePrefix {params.outdir_fprefix}/{params.run_id}."
        " --outSAMtype {params.output_SAM_type}"
        " --outSAMunmapped {params.output_SAM_unmapped}"
        " --outSAMattributes {params.output_SAM_attributes}"
        " {params.extra_args}"



#if Samtools_Merge_configs["do_merge"]:
#    rule Samtools_Merge:
#        input:
#            alignments = expand(rules.STAR_Align_Reads.output.alignment, sample_id=sample_ids)
#        resources:
#            walltime = Samtools_Merge_configs["walltime"],
#            nodes = Samtools_Merge_configs["nodes"],         
#            processors_per_node = Samtools_Merge_configs["processors_per_node"],
#            total_memory = Samtools_Merge_configs["total_memory"],  
#            logdir = "logs/",
#            job_id = (lambda wildcards: wildcards.sample_id)
#        shell:
#            "echo {input.alignments} > tomerge.txt"


##Mark Duplicates using Picard tools
