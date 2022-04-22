# File Name: RNA_Seq_Pipeline.snakefile
# Created By: ZW
# Created On: 2022-02-17
# Purpose: Manages the workflow for RNA-seq alignment and quantification for
# both single-end and paired-end sequencing


# library + module imports
# ----------------------------------------------------------------------------
import pandas as pd
from gardnersnake import ConfigurationHelper
from gardnersnake.misc import read_manifest_txt

# workflow setup
# ----------------------------------------------------------------------------

# build ConfigurationHelper and set the work dir
CH = ConfigurationHelper(cfg_dict=config)
analysis_id = CH.get_global_param("analysis_id")
workdir: CH.get_global_param("workdir", ispath=True, pathtype="dir", exists=True)


# workflow
#----------------

# Rule 0: Definig global outputs
# ----------------------------------------------------------------------------
rule all:
    input:
        star_cgi_rc = "star_create_genome_index.rc.out"


# Rule 1: Verify Genome Index or Create a new one.
# ----------------------------------------------------------------------------

# Collect info about where index should be located and what should be inside 
genome_index_path = CH.get_rule_param("STAR_CGI", "STAR_indices_path",
                        ispath=True, pathtype="dir", returnPath=True)
genome_index_manifest = read_manifest_txt(
                            CH.get_rule_param(rule="STAR_CGI",
                                param="index_manifest_file",
                                ispath=True, pathtype="file",
                                exists=True, returnPath=False
                            )
                        )

# check if the index already exists and verify its contents
if genome_index_path.exists():
    rule Verify_Index_Contents:
        input: index_dir = str(genome_index_path)
        params: manifest = genome_index_manifest
        output: star_cgi_rc = "star_create_genome_index.rc.out"
        resources: **CH.get_rule_resources("Verify_Index_Contents")
        shell:
            "check_directory --strict"
            " -o {output.star_cgi_rc}"
            " {params.manifest} {input.index_dir}"

# create the index if it doesn't already exist. verify contents afterwards
else:
    rule STAR_Create_Genome_Index:
        params:            
            nthreads = CH.get_rule_param(rule="STAR_CGI", param="nthreads"),
            sjdbOverhang = CH.get_rule_param(rule="STAR_CGI",param="sjdbOverhang"),
            manifest = genome_index_manifest,
            transcript_gtf_path = CH.get_rule_param(rule="STAR_CGI",
                                    param="transcript_gtf",
                                    ispath=True, pathtype="file",
                                    exists=True, returnPath=False
                                  ),
            fasta_path = CH.get_rule_param(rule="STAR_CGI",param="ref_fasta",
                            ispath=True, pathtype="file",
                            exists=True, returnPath=False
                         )
        resources: **CH.get_rule_resources("STAR_CGI") 
        output:
            genome_index_path = directory(str(genome_index_path)),
            star_cgi_rc = "star_create_genome_index.rc.out"
        shell:
            "mkdir -p {output.genome_index_path} && "
            "STAR --runThreadN {params.nthreads}"
            " --runMode genomeGenerate"
            " --genomeDir {output.genome_index_path}"
            " --genomeFastaFiles {params.fasta_path}"
            " --sjdbGTFfile {params.transcript_gtf_path}"
            " --sjdbOverhang {params.sjdbOverhang}"
            " && check_directory --strict"
            " -o {output.star_cgi_rc}"
            " {params.manifest}"
            " {output.genome_index_path}"


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



