#!python3

import os
import sys
import logging
import json
import copy
from MapTnSeq import RunMapTnSeq

def test():
    logging.basicConfig(level=logging.DEBUG)
    # Program should be run test.py config_fp.json fastq_fp genome_fna output_fp
    args = sys.argv
    config_fp = args[1]
    fastq_fp = args[2]
    genome_fna = args[3]
    output_fp = args[4]
    with open(config_fp, "r") as g:
        config_d = json.loads(g.read())

    program_d = config_d['values']
    program_d['fastq_fp'] = fastq_fp
    program_d['genome_fp'] = genome_fna
    program_d['output_fp'] = output_fp
    
    models_dir = os.path.join(os.getcwd(),"DATA/models")
    # Checking which model(s) are correct:
    succesful_models = FindBestModel(program_d, models_dir)

    # running program on model:
    program_d['modeltest'] = False
    program_d['maxReads'] = 10**8
    for mod in succesful_models:
        logging.info("Running 10**8 reads MapTnSeq on {}\n".format(mod))
        program_d['model_fp'] = mod
        RunMapTnSeq(program_d)

    


    return None


def FindBestModel(program_dict, models_dir):
    """
    Given the directory with all the models, run through each and do a 
    maxReads = 1000 test on each.
    models_dir is "./DATA/models"
    """

    # Getting all models files
    models_list = [os.path.join(models_dir, x) for x in os.listdir(models_dir)]
    # Ensuring max reads is 1000 
    program_dict['maxReads'] = 1000
    program_dict['modeltest'] = True
    succesful_models = []
    # Debug
    logging.debug(os.listdir(models_dir))
    # Run through all the models
    for model_fp in models_list:
        new_program_d = copy.deepcopy(program_dict)
        new_program_d['model_fp'] = model_fp
        returnCode = RunMapTnSeq(new_program_d)
        if isinstance(returnCode, int):
            if returnCode == 0:
                succesful_models.append(model_fp)
            logging.info("RETURN CODE: " + str(returnCode))

    logging.info("SUCCESFUL MODELS: " + "\n".join(succesful_models) + "\n")

    return succesful_models






def main():
    test()

    return None

if __name__ == "__main__":
    main()
