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
from pandas import DataFrame

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

# globals
# --------------------------
ALIGN_DIR = CH.get_parameters("STAR_Align_Reads")["outdir_fprefix"]
RUN_IDS = DM.get_shared_data({}, "run_name")

# workflow
#----------------

# rule 0: Defining global outputs
# ----------------------------------------------------------------------------
rule all: # static output unpacking should all go last
    input:
        expand(ALIGN_DIR + "{run_id}.Aligned.sortedByCoord.out.bam", run_id=RUN_IDS),
        expand(ALIGN_DIR + "{run_id}.Aligned.sortedByCoord.MarkDups.bam", run_id=RUN_IDS),
        expand(ALIGN_DIR + "{run_id}.MarkDups.metrics.txt", run_id=RUN_IDS),
        **DM.get_rule_data("STAR_Create_Genome_Index",["static_outputs"])


# rule 1: Verify Genome Index or Create a new one.
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


# rule 2: Aligning RNA-Seq reads using STAR
#-----------------------------------------------------------------------------
# scatter individual runs and align
rule STAR_Align_Reads:
    input:
        fastq1 = (lambda wildcards: DM.get_shared_data({'run_name': f"{wildcards.run_id}"},"fastq1")),
        fastq2 = (lambda wildcards: DM.get_shared_data({'run_name': f"{wildcards.run_id}"},"fastq2")),
        **DM.get_rule_data("STAR_Create_Genome_Index", ["static_outputs"])
    params:
        **CH.get_parameters("STAR_Align_Reads"),
        run_id = (lambda wildcards: wildcards.run_id),
        genome_index_dirn = CH.get_parameters("STAR_Create_Genome_Index")["genome_index_dir"]
    resources:
        ** CH.get_resources("STAR_Align_Reads", return_job_id=False),
        job_id = (lambda wildcards: wildcards.run_id)
    output: ALIGN_DIR + "{run_id}.Aligned.sortedByCoord.out.bam"
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

# merge output alignments to shared metadata file
align_table = DataFrame(zip(DM.get_shared_data({}, "run_name"),
    expand(ALIGN_DIR + "{run_id}.Aligned.sortedByCoord.out.bam", run_id=RUN_IDS)
    ), columns=["run_name", "alignment"])
DM.loj_shared_data(to_add=align_table, on=["run_name"], indicator=False)


# rule 3: MarkDuplicates in the alignments
# -----------------------------------------------------------------------------
# use gatk-picardtools to mark duplicates in alignments
rule GATK_Mark_Duplicates:
    input: (lambda wildcards: DM.get_shared_data({'run_name': f"{wildcards.run_id}"},"alignment"))
    params: **CH.get_parameters("GATK_Mark_Duplicates")
    output:
        metrics = (ALIGN_DIR + "{run_id}.MarkDups.metrics.txt"),
        align_md = (ALIGN_DIR + "{run_id}.Aligned.sortedByCoord.MarkDups.bam")
    resources:
        **CH.get_resources("GATK_Mark_Duplicates", return_job_id=False),
        job_id = (lambda wildcards: wildcards.run_id)
    shell:
        "gatk MarkDuplicates"
        " --java-options '-Xms{resources.total_memory}m -Xmx{resources.total_memory}m'"
        " -I {input} -M {output.metrics}"
        " -O {output.align_md} {params.extra_args}"

# add references to output MarkDups alignment files to shared metadata file
markdups_table = DataFrame(zip(DM.get_shared_data({}, "run_name"),
    expand(ALIGN_DIR + "{run_id}.Aligned.sortedByCoord.out.bam", run_id=RUN_IDS)),
    columns=["run_name", "alignment_markdups"])
DM.loj_shared_data(to_add=markdups_table, on=["run_name"], indicator=False)








# rule 3: Merge outputs into a single BAM file (optional)
#-----------------------------------------------------------------------------
# check if we want to create a merged output BAM and proceed accordingly
# do_merge = CH.get_parameters("Samtools_Merge_Alignments")["do_merge"]
# if do_merge:
#     rule Samtools_Merge_Alignments:
#         input: **DM.get_shared_data({},"alignment_markdups")
#         resources: **CH.get_resources("Samtools_Merge_Alignments", return_job_id=True)
#         output: **DM.get_rule_data("Samtools_Merge_Alignments", ["static_outputs"])
#         shell:
#             "samtools merge -o {output.merged_alignment_name} {input}"
#
#
# ##Mark Duplicates using Picard tools
