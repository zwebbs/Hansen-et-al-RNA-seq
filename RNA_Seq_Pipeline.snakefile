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
STAR_Align_configs = utils.extract_rule_params("STAR_Align", config)
Samtools_Merge_configs = utils.extract_rule_params("Samtools_Merge", config)
Picard_MarkDups_configs = utils.extract_rule_params("Picard_MarkDups", config)

# sample metadata processing
sample_metadata = pd.read_csv(config["sample_metadata"], sep=' ')
sample_ids = sample_metadata["sample"].tolist()
print(sample_metadata, sample_metadata.shape)

# workflow

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
            walltime = STAR_CGI_configs["walltime"],
            nodes = STAR_CGI_configs["nodes"],  
            processors_per_node = STAR_CGI_configs["processors_per_node"],
            total_memory = STAR_CGI_configs["total_memory"],
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


#Alignment using STAR
rule STAR_Align_Reads:
    input:
        genome_index = rules.STAR_Create_Genome_Indices.output.genome_index_dir,
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
        walltime = STAR_Align_configs["walltime"],
        nodes = STAR_Align_configs["nodes"],         
        processors_per_node = STAR_Align_configs["processors_per_node"],
        total_memory = STAR_Align_configs["total_memory"],  
        logdir = "logs/",
        job_id = (lambda wildcards: wildcards.sample_id)
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


if Samtools_Merge_configs["do_merge"]:
    rule Samtools_Merge:
        input:
            alignments = expand(rules.STAR_Align_Reads.output.alignment, sample_id=sample_ids)
        resources:
            walltime = Samtools_Merge_configs["walltime"],
            nodes = Samtools_Merge_configs["nodes"],         
            processors_per_node = Samtools_Merge_configs["processors_per_node"],
            total_memory = Samtools_Merge_configs["total_memory"],  
            logdir = "logs/",
            job_id = (lambda wildcards: wildcards.sample_id)
        shell:
            "echo {input.alignments} > tomerge.txt"


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



