# File Name: exceptions.py
# Created By: ZW
# Created On: 2022-02-24
# Purpose: defines user exceptions for snakemake workflow

# error class definitions

class UserError(Exception):
    ## base class for user-defined errors
    def __init__(self, message="Unspecified User Error"):
        self.message = message
        super().__init__(self.message)

