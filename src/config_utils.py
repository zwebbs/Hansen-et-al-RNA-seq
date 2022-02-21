# File Name: config_utils.py
# Created By: ZW
# Created On: 2022-02-21
# Purpose: Provides auxillary methods for readability and
# portability of snakemake file


# error class definitions

class UserError(Exception):
    ## base class for user-defined errors
    def __init__(self, message="Unspecified User Error"):
        self.message = message
        super().__init__(self.message)

# other class definitions


# method definitions

def extract_rule_params(rule_id, config_dict):
    # check that a valid parameters dictionary exists
    # and return it from original config dict
    rule_params = config_dict["rule_params"]
    if type(rule_params[rule_id]) is dict:
        return rule_params[rule_id]
    else:
        err_msg = f"Error. rule params could not be extracted "
        rule_msg = f"passed rule_id: ['{rule_id}'] "
        config_msg = f"check config dict: {config_dict}"
        raise UserError(message=err_msg+rule_msg+config_msg)



if __name__ == "__main__":
    pass
