#!python3

import os
import sys
import logging
import json
import random
from MapTnSeq import ParseModel



def generate_fake_inputs(inp_dict):
    """
    inp_dict must contain:
        model: (str) of model
        pastEnd: (str) of pastEnd
        num_tot: (int) number of total generated fake lines



    Needs to generate num_tot FASTQ reads.
    Random number generator will decide how many are usable,
        how many are not. (say num_tot = 100, usable = 60, unusable=40)
    For each real read, Also generate a fake quality.

    """





def test():
    args = sys.argv
    model_fp = args[1]
    num_tot = args[2]

    return None


def main():
    test()

    return None

if __name__ == "__main__":
    main()
