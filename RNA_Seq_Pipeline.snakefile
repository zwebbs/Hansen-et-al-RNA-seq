# File Name: RNA_Seq_Pipeline.snakefile
# Created By: ZW
# Created On: 2022-02-17
# Purpose: Manages the workflow for RNA-seq alignment and quantification for
# both single-end and paired-end sequencing

# library + module imports
import src.config_utils as utils


# configuration and setup
configfile: "config/RNA_Seq_PE_config_test.json"
STAR_CGI_configs = utils.extract_rule_params("STAR_CGI", config)
STAR_Align_configs = utils.extract_rule_params("STAR_align",config)


# workflow
rule all:
    input:
        directory(STAR_CGI_configs["STAR_indices_path"])


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
        shell:
            "STAR --runThreadN {params.nthreads}"
            " --runMode genomeGenerate"
            " --genomeDir {output.genome_index_dir}"
            " --genomeFastaFiles {params.fasta_path}"
            " --sjdbGTFfile {params.transcript_gtf_path}"
            " --sjdbOverhang {params.junction_overhang_limit}"


# Alignment using STAR
rule STAR_Align_Reads:
    input:
        genome_index_dir = STAR_CGI_configs["STAR_indices_path"]
        fastq_reads1 = utils.resolve_SEvPE_fastq(sequencing_metadata)
        fastq_reads2 = utils.resolve_SEvPE_fastq(sequencing_metadata)
    shell:
        "STAR --runThread {params.nthreads}"
        " --genomeDir {input.genome_index_dir} "
        ""
