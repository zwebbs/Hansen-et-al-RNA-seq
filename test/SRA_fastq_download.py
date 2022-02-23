# Fiile Name: SRA_fastq_download.py
# Created By: ZW
# Created On: 2021-02-170
# Purpose: Automate downloads of FASTQ files from NCBI SRA
# inspired by NCBI code documentation

import argparse
import subprocess


# Function Definitions
#---------------------

def get_SRA_ID_from_list(accession_list_file):
    # open a file and grab the SRA numbers from within
    # ONE SRA ID per line only.
    with open(accession_list_file, 'r') as fobj:
        sra_numbers = [line.strip() for line in fobj]
        return sra_numbers

def resolve_download_location(path, flag, alternative):
    if path is not None:
        return f"{flag} {path} "
    else: 
        return alternative
def resolve_extra_args(extra_args, alternative):
    if extra_args is not None:
        return f"{extra_args} "
    else:
        return alternative


def prefetch_fastq(sra_numbers, download_where=None, extra_args=None):
    # resolve flags for the commandline tool 'prefetch'
    outdir = resolve_download_location(download_where,"-O", "")
    prefetch_extra_args = resolve_extra_args(extra_args, "")
    # iterate through SRA_ids and run prefetch in detached subprocesses 
    for sra_id in sra_numbers:
        prefetch_cmd = f"prefetch {prefetch_extra_args} {outdir} {sra_id}"
        print (f"Currently Downloading: {sra_id}")
        print(f"The Command Was: {prefetch_cmd}")
        subprocess.call(prefetch_cmd, shell=True)
    print("Prefetch Command Complete. *(check output for errors)")
    return


## TODO currently throws an error if the fastq file already exists..
## solution might take some creativity because I'm not sure how consistently
## the files will be named here, or how to predict how many sets of files are expected.
def extract_fastq_from_sra(sra_numbers, input_dir=None, output_dir=None, extra_args=None):
    # loop through the SRA numbers and unpack the downloaded FASTQs from their
    # corresponding .sra files. also uses independent subprocess calls
    outdir_flag = resolve_download_location(output_dir,"--outdir", "")
    outdir = output_dir if output_dir is not None else ""
    fasterq_extra_args = resolve_extra_args(extra_args,"")
    for sra_id in sra_numbers:
        input_path = resolve_extra_args(input_dir, "").strip()
        sra_filepath = f"{input_path}{sra_id}/{sra_id}.sra" 
        fasterq_dump_cmd = f"fasterq-dump {outdir_flag} {fasterq_extra_args} {sra_filepath}" 
        gzip_fastq_cmd = f"gzip -v {output_dir}*.fastq"
        extract_cmd = fasterq_dump_cmd + " && " + gzip_fastq_cmd
        print(f"Currently Unpacking: {sra_id}")
        print(f"The Command Was: {extract_cmd}")
        subprocess.call(extract_cmd, shell=True)
    

# Main Execution
#----------------

if __name__ == "__main__":
    # setup commandline arguments for ease of use.
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefetch-extra-args", type=str)
    parser.add_argument("--prefetch-output-dir",type=str)
    parser.add_argument("--fasterq-dump-output-dir",type=str)
    parser.add_argument("--fasterq-dump-extra-args", type=str)
    parser.add_argument("sra_list", type=str)

    # parse arguments passed to the commandline
    args = parser.parse_args()
    print(args)
    # run the prefetch process on our samples
    sra_numbers = get_SRA_ID_from_list(args.sra_list)
    prefetch_fastq(sra_numbers=sra_numbers,
                    download_where=args.prefetch_output_dir,
                    extra_args=args.prefetch_extra_args)
    
    # extract fastq from SRA archive files and gzip
    extract_fastq_from_sra(sra_numbers=sra_numbers,
                            input_dir=args.prefetch_output_dir,
                            output_dir=args.fasterq_dump_output_dir,
                            extra_args=args.fasterq_dump_extra_args)




    exit(0)






