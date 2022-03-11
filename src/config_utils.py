# File Name: config_utils.py
# Created By: ZW
# Created On: 2022-02-21
# Purpose: Provides auxillary methods for readability and
# portability of snakemake file configuration

from exceptions import UserError
import pandas as pd

# method definitions

def extract_rule_params(rule_id, config_dict):
    # check that a valid parameters dictionary exists
    # and return it from original config dict
    rule_params = config_dict["rule_params"]
    if type(rule_params[rule_id]) is dict:
        return rule_params[rule_id]
    else:
        err_msg = "Error. rule params could not be extracted "
        rule_msg = f"passed rule_id: ['{rule_id}'] "
        config_msg = f"check config dict: {config_dict}"
        raise exceptions.UserError(message=err_msg+rule_msg+config_msg)


def STAR_manifest_from_sample_metadata(sample_metadata_tbl, outfile_path):
    # generate the STAR aligner file manifest from our sample metadata
    fastq_1, fastq_2, RG_info = [],[],[] # create empty lists for manifest
    for idx, dat in sample_metadata_tbl.iterrows():
        # generate fastq file lists automatically based on PE or SE
        sample_fastq1_list = dat["fastq1_pathlist"].split(',')
        fastq_1.extend(sample_fastq1_list)
        if dat["is_paired"] == 1:
            sample_fastq2_list = dat["fastq2_pathlist"].split(',')
            fastq_2.extend(sample_fastq2_list)
        else:
            fastq_2.extend(['-']*len(sample_fastq1_list))

        # generate sample ID list
        sample_id = dat["sample"]
        RG_info.extend([f"ID:{sample_id}"]*len(sample_fastq1_list))

    # generate and write final DataFrame
    fmtd_dat = [[f1,f2,rg] for (f1,f2,rg) in zip(fastq_1, fastq_2, RG_info)]
    fmtd_dat.to_csv(outfile_path, sep='\t',header=False, index=False)

if __name__ == "__main__":
    sample_metadata = pd.read_csv("test/sample_metadata_test.txt",sep="\t")
    STAR_manifest_from_sample_metadata(sample_metadata, "test/test_manifest.txt")
