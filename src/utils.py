# File Name: utils.py
# Created By: ZW
# Created On: 2022-02-21
# Purpose: Provides auxillary methods for readability and
# portability of snakemake workflow configurations

from .exceptions import UserError
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


def get_fastq_gz(sample_id, sample_metadata):
    # using the sample_id, grab the fastq data and configure 
    # the appropriate STAR command input for readFilesIn
    sample_indexed_meta = sample_metadata.set_index("sample")
    sample_dat = sample_indexed_meta.loc[sample_id]
    
    fastq1_pathlist = sample_dat["fastq1_pathlist"]
    fastq2_pathlist = sample_dat["fastq2_pathlist"]

    if sample_dat["is_paired"] == 1: 
        return f"{fastq1_pathlist} {fastq2_pathlist}"
    else:
        return f"{fastq1_pathlist}"

if __name__ == "__main__":
    sample_metadata = pd.read_csv("test/sample_metadata_test.txt", sep="\s+")
    print(get_fastq_gz("SRR17999770",sample_metadata))
