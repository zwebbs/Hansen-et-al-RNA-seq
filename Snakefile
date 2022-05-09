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





# Rule 1: Verify Genome Index or Create a new one.
# ----------------------------------------------------------------------------

# check if the index already exists and verify its contents
if not Path(DM.get_rule_data(rule_name="Verify_Index_Contents",
        key_list=["inputs", "genome_index_dir"])).exists():
    rule Verify_Index_Contents:
        input: **DM.get_rule_data("Verify_Index_Contents", ["inputs"])
        params: **CH.get_parameters("Verify_Index_Contents")
        output: **DM.get_rule_data("Verify_Index_Contents", ["outputs"])
        resources: **CH.get_resources("Verify_Index_Contents")
        shell:
            "check_directory --strict"
            " -o {output}"
            " {params.manifest} {input}"

# create the index if it doesn't already exist. verify contents afterwards
# else:  # create the index if it doesn't exist, verify contents afterwords
#     rule STAR_Create_Genome_Index:
#         input: **DM.get_rule_data("Verify_Index_Contents", ["input"])
#         params: **CH.get_parameters("Verify_Index_Contents")
#         output: **DM.get_rule_data("Verify_Index_Contents", ["output"])
#         resources: **CH.get_resources("Verify_Index_Contents")
#
#
#
#         shell:
#             "mkdir -p {output.genome_index_path} && "
#             "STAR --runThreadN {params.nthreads}"
#             " --runMode genomeGenerate"
#             " --genomeDir {output.genome_index_path}"
#             " --genomeFastaFiles {params.fasta_path}"
#             " --sjdbGTFfile {params.transcript_gtf_path}"
#             " --sjdbOverhang {params.sjdbOverhang}"
#             " && check_directory --strict"
#             " -o {output.star_cgi_rc}"
#             " {params.manifest}"
#             " {output.genome_index_path}"


# Rule 2: Aligning RNA-Seq reads using STAR
#-----------------------------------------------------------------------------


#TODO add handling for premade index paths
#Alignment using STAR
#rule STAR_Align_Reads:
#    input:
#        genome_index = rules.STAR_Create_Genome_Indices.output.genome_index_dir,
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
#rule Picard_MarkDups:
#    input:
#        alignment = STAR_Align_configs["outdir"] + "{sample_id}.Aligned.sortedByCoord.out.bam"
#    output:
#        alignment_dupsmarked = Picard_MarkDups_configs["outdir"] + "{sample_id}.MarkDups.sorted.bam",
#        markdups_metrics = Picard_MarkDups_configs["outdir"] + "{sample_id}.MarkDups.metrics.txt"
#    params:
#        extra_args = Picard_MarkDups_configs["extra_args"]
#    resources:
#        time = Picard_MarkDups_configs["cluster_time"],
#        nodes = Picard_MarkDups_configs["cluster_nodes"],
#        ntasks_per_node = Picard_MarkDups_configs["cluster_ntasks_per_node"],
#        cpus_per_task = Picard_MarkDups_configs["cluster_cpus_per_task"],
#        mem_per_cpu = Picard_MarkDups_configs["cluster_mem_per_cpu"],
#        logdir = "logs/",
#        job_id = lambda wildcards: wildcards.sample_id
#    shell:
#        "gatk MarkDuplicates"
#        " -I={input.alignment}"
#        " -O={output.alignment_dupsmarked}"
#        " -M={output.markdups_metrics}"
#        " -AS=true"
#        " {params.extra_args}"

#Filter Alignments by bitwise flags if specified 
#if run_filter_alignments:
#    rule Samtools_Filter_Alignments:
#        input:
#            pass
#



