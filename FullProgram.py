#!python3

import os
import sys
import logging
import json
from MapTnSeq import RunMapTnSeq
from DesignRandomPool import RunDesignRandomPool



def test():
    logging.basicConfig(level=logging.DEBUG)
    # Program should be run FullProgram.py map_cfg.json drp_cfg.json tmp_dir output_fp
    args = sys.argv
    map_cfg_fp = args[1] # map tn seq config file path (contains all)
    drp_cfg_fp = args[2] # design random pool config file path (contains all)
    tmp_dir = args[3] # tmp dir to do work
    pool_output_fp = args[4] # the pool file

    # Checking inputs exist
    for x in [map_cfg_fp, drp_cfg_fp]:
        if not os.path.exists(x):
            raise Exception(x + " does not exist, cannot run without configs")
    for x in [tmp_dir, os.path.dirname(pool_output_fp)]:
        if not os.path.isdir(x):
            raise Exception("Input directory " + x + " doesn't exist.")
    if os.path.exists(pool_output_fp):
        resp = input("Pool file output file already exists. Overwrite? y/n")
        if resp == "n":
            sys.exit(0)


    # Loading MapTNSEQ Config Dictionary from JSON to python
    with open(map_cfg_fp, "r") as g:
        map_cfg = json.loads(g.read())


    # Find the model:
    models_dir = os.path.join(os.getcwd(), "DATA/models")
    if not os.path.isdir(models_dir):
        raise Exception("Models Directory doesn't exist: " + models_dir)

    good_models_list = FindWorkingModel(map_cfg, models_dir)

    if len(good_models_list) > 1:
        logging.warning("Found more than one working model file: " \
                + "Using the first of these: \n" + "\n".join(good_models_list))
    elif len(good_models_list) == 0:
        raise Exception("No good models found in models dir: " \
                + "\n".join(os.listdir(models_dir)))
    else:
        logging.info("Found a single good model to use: " + good_models_list[0])

    model_fp = good_models_list[0]

    map_cfg["model_fp"] = model_fp

    num_mts_runs = len(map_cfg['fastq_fp_list'])
    if  num_mts_runs > 100:
        raise Exception("Cannot run Map Tnseq on more than 100 files " \
                + "\n Currently: " + str(num_mts_runs))


    MapTS_Output_fps = []
    for i in range(len(map_cfg['fastq_fp_list'])):
        map_cfg["fastq_fp"] = map_cfg['fastq_fp_list'][i]

        # Here we write the output filepath of the program 
        if i < 10:
            suffix = "0"
        else:
            suffix = ""

        # (current MapTnSeq Table File)
        cMTS_output_fp = os.path.join(tmp_dir, "MTS_Table_" + suffix + str(i))
        map_cfg['output_fp'] = cMTS_output_fp
        MapTS_Output_fps.append(cMTS_output_fp)

        RunMapTnSeq(map_cfg)

    # LOADING DESIGN RANDOM POOL INPUT
    with open(drp_cfg_fp, "r") as g:
        drp_cfg = json.loads(g.read())

    # Adding map tn seq table files:
    drp_cfg["map_tnseq_table_fps"] = MapTS_Output_fps
    drp_cfg["output_fp"] =  pool_output_fp

    # Running Design Random Pool
    RunDesignRandomPool(drp_cfg)


    return None


def FindWorkingModel(program_dict, models_dir):
    """
    Given the directory with all the models, run through each and do a 
    maxReads = 1000 test on each.
    models_dir is "./DATA/models"
    We need a single fastq file out of the list to do the test run.
    We choose the first one in the list fastq_fp_list
    """

    # Getting all models files
    models_list = [os.path.join(models_dir, x) for x in os.listdir(models_dir)]

    # Getting the first fastq file set
    program_dict['fastq_fp'] = program_dict['fastq_fp_list'][0]

    # Ensuring max reads is 1000 
    program_dict['maxReads'] = 1000
    program_dict['modeltest'] = True
    succesful_models = []

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
