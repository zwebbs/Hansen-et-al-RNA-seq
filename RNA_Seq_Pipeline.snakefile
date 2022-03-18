# File Name: RNA_Seq_Pipeline.snakefile
# Created By: ZW
# Created On: 2022-02-17
# Purpose: Manages the workflow for RNA-seq alignment and quantification for
# both single-end and paired-end sequencing

# library + module imports
import pandas as pd
import src.utils as utils

# configuration and setup
# configfile: is not hardcoded and must be supplied through the CLI
analysis_id = config["analysis_id"]
STAR_CGI_configs = utils.extract_rule_params("STAR_CGI", config)
STAR_Align_configs = utils.extract_rule_params("STAR_align",config)

# sample metadata processing
sample_metadata = pd.read_csv(config["sample_metadata"], sep='\s+')
sample_ids = sample_metadata["sample"].tolist()


# workflow
rule all:
    input:
        directory(STAR_CGI_configs["STAR_indices_path"]),
        expand(STAR_Align_configs["outdir"] + "{sample_id}.Aligned.sortedByCoord.out.bam", sample_id=sample_ids)


# Create Genome Index using STAR if not already created
if STAR_CGI_configs["premade_index_path"] is None:
    print("Generating Index for STAR...")
    rule STAR_Create_Genome_Indices:
        output:
            genome_index_dir = directory(STAR_CGI_configs["STAR_indices_path"])
        params:
            nthreads = STAR_CGI_configs["nthreads"],
            fasta_path = STAR_CGI_configs["ref_fasta"],
            transcript_gtf_path = STAR_CGI_configs["transcript_gtf"],
            junction_overhang_limit = STAR_CGI_configs["sjdbOverhang"]
        resources:
            time = STAR_CGI_configs["cluster_time"],
            nodes = STAR_CGI_configs["cluster_nodes"],  
            ntasks_per_node = STAR_CGI_configs["cluster_ntasks_per_node"],
            cpus_per_task = STAR_CGI_configs["cluster_cpus_per_task"],
            mem_per_cpu = STAR_CGI_configs["cluster_mem_per_cpu"],
            logdir = "logs/",
            job_id = ""
        shell:
            "mkdir {output.genome_index_dir} && "
            "STAR --runThreadN {params.nthreads}"
            " --runMode genomeGenerate"
            " --genomeDir {output.genome_index_dir}"
            " --genomeFastaFiles {params.fasta_path}"
            " --sjdbGTFfile {params.transcript_gtf_path}"
            " --sjdbOverhang {params.junction_overhang_limit}"


# Alignment using STAR
rule STAR_Align_Reads:
    input:
        genome_index_dir = STAR_CGI_configs["STAR_indices_path"],
    params:
        sample_id = (lambda wildcards: wildcards.sample_id),
        sample_fastq_gz = (lambda wildcards: utils.get_fastq_gz(wildcards.sample_id,sample_metadata)),
        nthreads = STAR_Align_configs["nthreads"],
        outdir_fprefix = STAR_Align_configs["outdir"],
        output_SAM_type = STAR_Align_configs["outSAMtype"],
        output_SAM_unmapped = STAR_Align_configs["outSAMunmapped"],
        output_SAM_attributes = STAR_Align_configs["outSAMattributes"],
        extra_args = STAR_Align_configs["extra_args"]
    output:
        alignment = STAR_Align_configs["outdir"] + "{sample_id}.Aligned.sortedByCoord.out.bam"
    resources:
        time = STAR_Align_configs["cluster_time"],
        nodes = STAR_Align_configs["cluster_nodes"],         
        ntasks_per_node = STAR_Align_configs["cluster_ntasks_per_node"],
        cpus_per_task = STAR_Align_configs["cluster_cpus_per_task"],  
        mem_per_cpu = STAR_Align_configs["cluster_mem_per_cpu"],
        logdir = "logs/",
        job_id = lambda wildcards: wildcards.sample_id
    shell:
        "STAR --runThreadN {params.nthreads}"
        " --runMode alignReads"
        " --quantMode GeneCounts"
        " --genomeDir {input.genome_index_dir}/"
        " --readFilesIn {params.sample_fastq_gz}"
        " --readFilesCommand gunzip -c"
        " --outFileNamePrefix {params.outdir_fprefix}/{params.sample_id}."
        " --outSAMtype {params.output_SAM_type}"
        " --outSAMunmapped {params.output_SAM_unmapped}"
        " --outSAMattributes {params.output_SAM_attributes}"
        " {params.extra_args}"


