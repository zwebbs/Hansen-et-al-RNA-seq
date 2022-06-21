# File Name: get-cisregulatory-genes.sh
# Created By: ZW
# Created On: 2022-05-19
# Purpose: given SNP positions, extract genes from NCBI gtf which are within
# the passed genome window, and return a file with SNP the SNP and corresponding
# genes

# check passed commandline arguments
## of which there are three.
###   -s <snp info BED6 file .bed> [required]
###   -g <NCBI transcript/gene models .gtf> [required]
###   -w (INT genome window -number of base pairs upstream and downstream
###       to define the cis-regulatory region of interest) [required]

while getopts ":s:g:w" 'opt';
do
    case ${opt} in
        s) 
            snpinfo_file="${OPTARG}"
            ;; # capture the SNPinfo file (3 columns) BED6 format

        g) 
            gtf_file="${OPTARG}"
            ;; # capture the gtf file containing the gene models

        w)
            window_nt=${OPTARG}
            ;; # capture the genome window size in bp that defines the
               #cis-regulatory region
        *)
            echo "INVALID_OPTION -- ${OPTARG}"
            exit 1
            ;; # capture all other (invalid) inputs
    esac
done

# intersect the window around each SNP with a list of genes
while IFS=, read chrom start end rsID other; do
    printf '"%s"\t"%s"\t"%s"\t"%s"\n' "${chrom}" "${start}" "${end}" "${rsID}"
done <<<$(cat ${snpinfo_file})




#bedtools window -w 1000000 -a snpinfo.bed -b hg38.p13.NCBI_RefSeq.gtf
