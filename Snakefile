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


# workflow
#----------------

# Rule 0: Defining global outputs
# ----------------------------------------------------------------------------

rule all:
    input:
        **DM.get_rule_data("Verify_Index_Contents", ["outputs"])


# Rule 1: Verify Genome Index or Create a new one.
# ----------------------------------------------------------------------------

# check if the index already exists and verify its contents
genome_index_dir = CH.get_parameters("Verify_Index_Contents")["genome_index_dir"]
if Path(genome_index_dir).exists():
    rule Verify_Index_Contents:
        params: **CH.get_parameters("Verify_Index_Contents")
        output: **DM.get_rule_data("Verify_Index_Contents", ["outputs"])
        resources: **CH.get_resources("Verify_Index_Contents")
        shell:
            "check_directory --strict"
            " -o {output}"
            " {params.genome_index_manifest} {params.genome_index_dir}"

else: # create the index if it doesn't exist, verify contents afterwords
    rule STAR_Create_Genome_Index:
        input: **DM.get_rule_data("STAR_Create_Genome_Index", ["inputs"])
        params:
            **CH.get_parameters("Verify_Index_Contents"),
            **CH.get_parameters("STAR_Create_Genome_Index")
        output: **DM.get_rule_data("Verify_Index_Contents", ["outputs"])
        resources: **CH.get_resources("STAR_Create_Genome_Index")
        shell:
            "mkdir -p {params.genome_index_dir} && "
            "STAR --runThreadN {params.nthreads}"
            " --runMode genomeGenerate"
            " --genomeDir {params.genome_index_dir}"
            " --genomeFastaFiles {input.genome_fasta_path}"   
            " --sjdbGTFfile {input.transcript_gtf_path}"
            " --sjdbOverhang {params.sjdbOverhang}"
            " && check_directory --strict"
            " -o {output}"
            " {params.genome_index_manifest}"
            " {params.genome_index_dir}"


# Rule 2: Aligning RNA-Seq reads using STAR
#-----------------------------------------------------------------------------

# split alignment by sequencing runs
rule STAR_Align_Reads:
    input:
        **DM.get_rule_data("Verify_Index_Contents", ["outputs"]),
        fastqs=DM.get_from_shared_data(["libraries", "runs", "run_name"],"{run_name}",["fastq1","fastq2"])
#    params:
#        sample_id = (lambda wildcards: wildcards.sample_id),
#        sample_fastq_gz = (lambda wildcards: utils.get_fastq_gz(wildcards.sample_id,sample_metadata)),
#        nthreads = STAR_Align_configs["nthreads"],
#        outdir_fprefix = STAR_Align_configs["outdir"],
#        output_SAM_type = STAR_Align_configs["outSAMtype"],
#        output_SAM_unmapped = STAR_Align_configs["outSAMunmapped"],
#        output_SAM_attributes = STAR_Align_configs["outSAMattributes"],
#        extra_args = STAR_Align_configs["extra_args"]
#    output:
#        alignment = STAR_Align_configs["outdir"] + "{sample_id}.Aligned.sortedByCoord.out.bam"
#    resources:
#        walltime = STAR_Align_configs["walltime"],
#        nodes = STAR_Align_configs["nodes"],         
#        processors_per_node = STAR_Align_configs["processors_per_node"],
#        total_memory = STAR_Align_configs["total_memory"],  
#        logdir = "logs/",
#        job_id = (lambda wildcards: wildcards.sample_id)
#    shell:
#        "STAR --runThreadN {params.nthreads}"
#        " --runMode alignReads"
#        " --quantMode GeneCounts"
#        " --genomeDir {input.genome_index_dir}/"
#        " --readFilesIn {params.sample_fastq_gz}"
#        " --readFilesCommand gunzip -c"
#        " --outFileNamePrefix {params.outdir_fprefix}/{params.sample_id}."
#        " --outSAMtype {params.output_SAM_type}"
#        " --outSAMunmapped {params.output_SAM_unmapped}"
#        " --outSAMattributes {params.output_SAM_attributes}"
#        " {params.extra_args}"



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
