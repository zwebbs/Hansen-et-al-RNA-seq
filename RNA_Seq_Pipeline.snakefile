# File Name: RNA_Seq_Pipeline.snakefile
# Created By: ZW
# Created On: 2022-02-17
# Purpose: Manages the workflow for RNA-seq alignment and quantification for
# both single-end and paired-end sequencing

# library + module imports
import pandas as pd
import src.config_utils as conf_utils

# configuration and setup
# configfile: is not hardcoded and must be supplied through the CLI
analysis_id = config["analysis_id"]
STAR_CGI_configs = conf_utils.extract_rule_params("STAR_CGI", config)
STAR_Align_configs = conf_utils.extract_rule_params("STAR_align",config)

# sample metadata processing
sample_metadata = pd.read_csv(config["sample_metadata"])
STAR_manifest_filename = f"{analysis_id}.manifest"
STAR_manifest_from_sample_metadata(sample_metadata, STAR_manifest_filename)

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
        file_manifest = STAR_manifest_filename
    params:
        nthreads = STAR_Align_configs["nthreads"],
        outdir_fprefix = STAR_Align_configs["outdir"] + analysis_id,
        output_SAM_type = STAR_Align_configs["outSAMtype"],
        output_SAM_unmapped = STAR_Align_configs["outSAMunmapped"],
        output_SAM_attributes = STAR_Align_configs["outSAMattributes"],
        extra_args = STAR_Align_configs["extra_args"]
    shell:
        "STAR --runThread {params.nthreads}"
        " --genomeDir {input.genome_index_dir}"
        " --readFilesManifest {input.file_manifest}"
        " --outFileNamePrefix {params.outdir_fprefix}"
        " --outSAMtype {params.output_SAM_type}"
        " --outSAMunmapped {params.output_SAM_unmapped}"
        " --outSAMattributes {params.output_SAM_attributes}"
        " {params.extra_args}"
