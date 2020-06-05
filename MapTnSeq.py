#!python3

import os
import sys
import json
import logging
import argparse
import re
import subprocess
import copy

"""
Notes:
    Program called through "RunMapTnSeq"
    Overall scratch dir will be scratch/MapTStmp
    temporary fna file will be called "temp.fna"
    trunc file? will be called "trunc.txt"
    FASTQ file will be unzipped.


    What is Unmapped file
    What is trunc file?
    

    input_args should be a dict with the following:
        debug: (bool)
        modeltest: (bool)
        maxReads: (int) or None
        minQuality: (int)
        flanking: (int)
        wobbleAllowed: (int)
        minIdentity: (int)
        minScore: (int)
        delta: (int)
        tileSize: (int)
        stepSize: (int)
        nMapped: (int)
        nMapUnique: (int)
        nPastEndIgnored: (int)
        nPastEndTrumps:  (int)
        tmp_dir: (str) Path to working directory
        blatcmd: (str) BLAT command
        unmapped_fp: (str) Path to write unmapped file
        trunc_fp: (str) Path to write TRUNC file
        genome_fp: (str) Path to genome FASTA file
        model_fp: (str) filepath to model 
        fastq_fp: (str) Path to fastq file
                Note, fastq file should be unzipped. If it contains 
                paired-end reads, the second read will be ignored.
        output_fp: (str) Path to output file

    Things to check:
    What does tmpFNA look like?
    What does input FASTQ file look like?
    What does model look like?
    Do you get what you expect?

    what is BLAT8? What does it do?



"""


def RunMapTnSeq(input_args):
    """
    This function runs the entire process of MapTnSeq
        The variable "output_fp" in input_args is where 
        the output of the program is written to.
        The output file is tab-delimited and contains a list
        of the usable reads, where, by column, the information
        is presented: read name, barcode, which scaffold, pos (bp?),
        strand, unique or not flag, beginning and end of hit to genome in 
        read after trimming the transposon seq, the bit score, and
        the % identity (10 overall fields)
    """
    
    # We check that inputs are correct and add "model" and "pastEnd":
    parsed_vars = ParseInputs(input_args)

    
    # generates the dict "mapnames", parses input fastq file and finds barcodes.
    ret_d = find_barcodes_and_end_of_transposon(parsed_vars)
    # It is possible nothing worked, and thus ret_d = 1
    if isinstance(ret_d, int) and ret_d == -1:
        logging.warning("Stopping current run on:\n"
        "Model {}\n".format(parsed_vars['model_fp']) + \
        "FASTQ {}\n".format(parsed_vars['fastq_fp']) + \
        "Genome {}".format(parsed_vars['genome_fp']))
        return -1
    else:
        parsed_vars.update(ret_d)
        if parsed_vars['modeltest']:
            logging.info("Model: {} passed".format(parsed_vars['model_fp']))
            return 0

    if parsed_vars['pastEnd'] is not None and parsed_vars['pastEnd'] != '':
        # generates the dict "hitsPastEnd", runs BLAT8 
        hitsPastEnd = pastEndBLAT8(parsed_vars)
    else:
        hitsPastEnd = {}

    parsed_vars['hitsPastEnd'] = hitsPastEnd
  

    # Creates Output File Handle, to be closed later.
    OUTPUT, HG_d = Map_to_genome(parsed_vars)
    parsed_vars['OUTPUT'] = OUTPUT
    parsed_vars.update(HG_d)


    # Writes the file "unmapped_fp" in parsed_vars
    # also, unlinks the file tmpFNA_fp if debug isn't True
    nUnmapped = WriteUnMappedFile(parsed_vars)
    parsed_vars['nUnmapped'] = nUnmapped

    # We add to "OUTPUT" file handle all the pastEnd values
    PrintOutHitsPastEnd(parsed_vars)

    parsed_vars['OUTPUT'].close()

    # Generate text report
    parsed_vars['text_report'] = CreateTextReport(
                    parsed_vars['nLong'], 
                    parsed_vars['nReads'],
                    parsed_vars 
                    )

    logging.info(parsed_vars['text_report'])

    logging.info("Finished Running MapTnSeq. Return dict empty.")

    return_dict = {
            }

    return return_dict



def ParseInputs(input_args):
    """
    input_args should be a dict with the following:
        debug: (bool)
        modeltest: (bool)
        maxReads: (int or None)
        minQuality: (int)
        flanking: (int)
        wobbleAllowed: (int)
        minIdentity: (int)
        minScore: (int)
        delta: (int)
        nMapped: (int)
        nMapUnique: (int)
        nPastEndIgnored: (int)
        nPastEndTrumps:  (int)
        tileSize: (int)
        stepSize: (int)
        tmp_dir: (str) Path to working directory
        blatcmd: (str) BLAT command
        unmapped_fp: (str) Path to write unmapped file
        tmpFNA_fp: (str) Path to Temp FNA file
        trunc_fp: (str) Path to write TRUNC file
        genome_fp: (str) Path to genome FASTA file
        model_fp: (str) filepath to model 
        fastq_fp: (str) Path to fastq file
                Note, fastq file should be unzipped. If it contains 
                paired-end reads, the second read will be ignored.
        output_fp: (str) Path to output file
    """

    # Checking debug and modeltest
    for v in ["debug", "modeltest"]:
        if v not in input_args or not isinstance(input_args[v], bool):
            raise Exception(v + "must be included in inputs and must be bool.")

    # Checking ints
    for v in ["maxReads", "minQuality", "flanking",
            "wobbleAllowed", "minIdentity", "minScore", "delta",
            "tileSize", "stepSize", "nMapped", "nMapUnique",
            "nPastEndIgnored", "nPastEndTrumps"]:
        if v not in input_args or not isinstance(input_args[v], int):
            if v == "maxReads" and input_args[v] is None:
                continue
            else:
                raise Exception(v + " must be in inputs and must be int.")

    # Checking strs
    for v in ["tmp_dir", "unmapped_fp", "trunc_fp", "genome_fp",
            "model_fp", "fastq_fp", "output_fp", "blatcmd", "tmpFNA_fp"]:
        if v not in input_args or not isinstance(input_args[v], str):
            raise Exception(v + " must be in inputs and must be str.")


    # Parsing Model File into "model" and "pastEnd":
    # Model
    model, pastEnd = ParseModel(input_args['model_fp'])

    input_args['model'] = model
    input_args['pastEnd'] = pastEnd

    return input_args


def ParseModel(model_fp):
    """
    Inputs:
        model_fp: (str) is path to model file, should be 1 or 2 lines
    Outputs:
        model: (str) the string of the first line of the model file.
        pastEnd: (str) string of the second line of the model file. Could be ''
    """
    with open(model_fp, "r") as g:
        model = g.readline()
        if not re.match(r'^n*[ACGT]+N+[ACGT]+$', model):
            raise Exception("Invalid model: " + model + "\n for file " + \
                    input_args['model_fp'])
        pastEnd = g.readline()
        if pastEnd != '':
            if not re.match(r'^[ACGT]+$', pastEnd):
                raise Exception("Invalid past-end sequence: " + pastEnd + \
                        "\n for file " + input_args['model_fp'])

    return [model, pastEnd]




def CreateTextReport(nLong, nReads, inp_dict):
    """
    We create a string with info about program:

    nReads (int)
    nLong (int)
    inp_dict: (dict) Contains:
        nameToBarcode (dict)
        nTryToMap (int) 
        nMapped (int)
        nMapUnique (int)
        hitsPastEnd (dict)
        nPastEndIgnored (int)
        nPastEndTrumps (int)
    """
    text_str = "Reads processed (nReads): " + str(nReads) + "\n" \
              + "Long-enough (nLong): " + str(nLong) + "\n" \
              + "Barcodes found (nameToBarcode): " + str(len(inp_dict['nameToBarcode'].keys())) + "\n" \
              + "Mapping attempted for (nTryToMap): " + str(inp_dict['nTryToMap']) + "\n" \
              + "Mapped (nMapped): " + str(inp_dict['nMapped']) + "\n" \
              + "Uniquely (nMapUnique): " + str(inp_dict['nMapUnique']) + "\n" \
              + "Hits past end of transposon (hitsPastEnd - nPastEndIgnored): " \
              + str(len(inp_dict['hitsPastEnd'].keys()) \
                        - inp_dict['nPastEndIgnored']) + "\n" \
                        + "Weak/ambiguous (nPastEndIgnored): " + str(inp_dict['nPastEndIgnored']) + "\n" \
              + "Trumped hit to genome (nPastEndTrumps):" + str(inp_dict['nPastEndTrumps']) + " times\n"
    if inp_dict['nReads'] > 0:
        text_str += "Proportions: \n" \
                  + "Long-enough (nLong/nReads): {:.3f}\n".format(nLong/nReads) \
                  + "Barcode (nameToBarcode/nReads): {:3f}\n".format(len(inp_dict['nameToBarcode'].keys())/nReads) \
                  + "Attempted (nTryToMap/nReads): {:3f}\n".format(inp_dict['nTryToMap']/nReads) \
                  + "Mapped (nMapped/nReads): {:3f}\n".format(inp_dict['nMapped']/nReads) \
                  + "Past-end (hits past end/ nReads): {:3f}\n".format(
                    (len(inp_dict['hitsPastEnd'].keys())-inp_dict['nPastEndIgnored'])/nReads)

    return text_str









def PrintOutHitsPastEnd(inp_dict):
    """
    This function writes "hits past end" to OUTPUT file handle

    Inputs:
        hitsPastEnd (dict)
        nameToBarcode (dict)
        OUTPUT: file handle to output

    """
    for k in inp_dict['hitsPastEnd']:
        read, score = k, inp_dict['hitsPastEnd'][k]
        if not score > 0:
            continue
        OUTPUT.write("\t".join([
            read, 
            str(inp_dict['nameToBarcode'][read]),
            "pastEnd",
            "",
            "",
            "",
            "",
            "",
            "",
            ""
            ]))






def WriteUnMappedFile(inp_dict):
    """
    Inputs:

    inp_dict: (dict) must contain:
        unmapped_fp: (str)
        tmpFNA_fp: (str)
        mapnames: (dict)
        nameToBarcode: (dict)
        debug: (bool)

    vars created:
        nUnmapped (int) number of lines in unMapped file

    Writing to unmapped file at "unmapped_fp", we check each barcode in 
        tmpFNA_fp and see if it's in the dict "mapnames"
    """

    nUnmapped = 0
    tmpFNA_fp = inp_dict['tmpFNA_fp']
    UNMAPPED = open(inp_dict['unmapped_fp'], "w")
    FNA = open(tmpFNA_fp, "r")
    line_num = 1
    header = FNA.readline()
    while header != "":
        header = header.rstrip()
        if not header[0] == ">":
            raise Exception("Cannot parse {} at line {} in {}".format(header, 
                line_num, tmpFNA_fp))
        name = header[1:]
        seq = FNA.readline()
        line_num += 1
        seq = seq.rstrip()
        if not re.match(r'[ACGTN]+$', seq):
            raise Exception("Cannot parse {} in {} at line {}".format(
                            seq,tmpFNA_fp, line_num))
        if name in inp_dict['mapnames'] and inp_dict['mapnames'][name] == 0:
            UNMAPPED.write('>{} {}\n{}\n'.format(
                inp_dict['nameToBarcode'][name],
                name,
                seq
                ))
            nUnmapped += 1
        header = FNA.readline()
        line_num += 1
    

    FNA.close()
    UNMAPPED.close()
    logging.info("Wrote {} unmapped reads to {} in fasta format\n".format(
        nUnmapped,
        inp_dict['unmapped_fp']))

    if not inp_dict['debug']:
        os.unlink(tmpFNA_fp)

    return nUnmapped







def Map_to_genome(inp_dict):
    """
    inp_dict:
        debug: (bool)
        unmapped_fp: (str)
        hitsPastEnd: (dict)
        mapnames: (dict)
        tmpFNA_fp: (str)
        genome_fp: (str) One of the inputs to the program (as genome) should
            be a FASTA file
        blatcmd: (str) (Location of blat8 in OS)
        tmp_dir: (str) path to working dir
        blatcmd: (str) Path to BLAT command, i.e. Bin/blat
        minScore: (int) minimum score for mapping to genome or past-end
        minIdentity: (int) minimum %identity for mappinf to genome or past-end
        tileSize: (int) size of an alignment tile
        stepSize: (int) distance between the starting bases of alignment tiles
            (will overlap if stepSize<tileSize)
        output_fp: (str) path to output_file

        The below are used for Handle Genome BLAT

        nPastEndTrumps: (int)  hit to past-end (close to) as good as hit to genome
        nPastEndIgnored: (int) weak hit to past-end ignored 
        nMapUnique: (int)
        nMapped: (int)
        nameToBarcode: (dict)
        delta: (int) minimum difference in score for considering a hit unique.
           
    Returns:
        OUTPUT file handle for output file

    """
    b8_d = inp_dict
    b8_d['queriesFile'] = inp_dict['tmpFNA_fp']
    b8_d['dbFile'] = inp_dict['genome_fp']
    blat8_fp = os.path.join(inp_dict['tmp_dir'], 'blat8_genome')
    b8_d['blat8_fp'] = blat8_fp
    RunBLAT8(b8_d)

    logging.info("Parsing " + blat8_fp)

    b = open(blat8_fp, "r")
    # OUTPUT is file handle for output
    OUTPUT = open(inp_dict['output_fp'], "w")

    HG_d = {
        "nPastEndTrumps": inp_dict[ "nPastEndTrumps"],
        "nPastEndIgnored": inp_dict[ "nPastEndIgnored"],
        "nMapUnique": inp_dict[ "nMapUnique"],
        "nMapped": inp_dict[ "nMapped"],
        "nameToBarcode": inp_dict[ "nameToBarcode"],
        "delta": inp_dict[ "delta"],
        "OUTPUT": OUTPUT
        }


    lines = []
    c_line = b.readline()

    while c_line != "":
        # Removing newline \n
        c_line = c_line.rstrip()
        F = c_line.split('\t')
        query, subject, identity, l, mm, gs, qB, qE, sB, sE, evl, score = F

        inp_dict['mapnames'][query] = 1

        if len(lines) == 0 or query == lines[0][0]:
            lines.append(F)
        else:
            HandleGenomeBLAT(lines, inp_dict['hitsPastEnd'], HG_d, 
                    inp_dict['debug'])
            lines = [F]
        c_line = b.readline()

    HandleGenomeBLAT(lines, inp_dict['hitsPastEnd'], HG_d, inp_dict['debug'])

    b.close()

    if not inp_dict['debug']:
        os.unlink(blat8_fp)

    # We return file handle, and HG_d
    return [OUTPUT, HG_d]
           










    





# This only runs if there is a "past End" part to the model (a second line)
def pastEndBLAT8(inp_dict):
    """
    Inputs:

    inp_dict: (dict) must contain
        endFNA_fp: (str) filepath to End FNA
        tmpFNA_fp: (str) filepath to temp FNA
        pastEnd: (str) an optional part of the end of the model
            that contains the sequence "past the end" of the 
            transposon that might arise from residual intact plasmid.
        blatcmd: (str) BLAT command
        tmp_dir: (str) path to working directory
        minScore: (int) minimum score for mapping to genome or past-end
        minIdentity: (int) minimum %identity for mappinf to genome or past-end
        tileSize: (int) size of an alignment tile
        stepSize: (int) distance between the starting bases of alignment tiles
            (will overlap if stepSize<tileSize)
        debug: (bool)
    Optional:
        mapnames: (dict)

    Outputs:
        hitsPastEnd: (dict) read (str) to score (float) of match to past-end sequence

    """

    # read to score of match to past-end sequence
    hitsPastEnd = {}

    with open(inp_dict['endFNA_fp'], "w") as g:
        g.write(">pastend\n{}\n".format(inp_dict['pastEnd']))

    #Making the input to RunBLAT8:
    blat8_inp_d = inp_dict
    blat8_inp_d['queriesFile'] = inp_dict['tmpFNA_fp']
    blat8_inp_d['dbFile'] = inp_dict['endFNA_fp']
    blat8_fp = os.path.join(inp_dict['tmp_dir'], 'blat8_pastEnd')
    blat8_inp_d['blat8_fp'] = blat8_fp 
    RunBLAT8(blat8_inp_d)

    if inp_dict['debug']:
        logging.info("Parsing past-end hits to {}\n".format(blat8_fp))

    b = open(blat8_fp, "r")
    crnt_line =  b.readline()
    while crnt_line != '':
        crnt_line = crnt_line.rstrip()
        F = crnt_line.split('\t')
        query, subject, identity, l, mm, gs, qB, qE, sB, sE, evl, score = F
        if not (query in hitsPastEnd and hitsPastEnd[query] > score):
            hitsPastEnd[query] = score
        if 'mapnames' in inp_dict:
            inp_dict['mapnames'][query] = 1
        hitsPastEnd[query] = float(scr)
        crnt_line = b.readline()

    b.close()

    #In case of debug, we keep the blat8 file for checking later.
    if not inp_dict['debug']:
        os.unlink(blat8_fp)

    return hitsPastEnd 







# This function finds barcodes and end of transposon and writes remaining 
# sequence to TMP FNA
# FASTQ must already be unzipped
def find_barcodes_and_end_of_transposon(inp_dict):
    """
    This function does the following:
        Reads FASTQ input file.


    Inputs:
    inp_dict must contain:
        fastq_fp (str) Read from
        trunc_fp (str) Write to
        tmpFNA_fp (str) Write to
        unmapped_fp (str or None) Write to (fp)
        maxReads (int or None) (maxReads) 
        model (str) Model string
        nameToBarcode: (dict)
        minScore (int) The minimum amount of sequence beyond model needed.
        wobbleAllowed (int) uncertainty in location of barcode or end of transposon,
            on either side of expectation
        minQuality (int) every nucleotide in a barcode must be at least this quality
        flanking (int) number of nucleotides on each side that must match
        debug (bool)

    Returns:
        ret_d: (dict) containing
            nLong
            nTryToMap
            nReads
            mapnames

    """

    limit = inp_dict['maxReads']
    nameToBarcode = inp_dict['nameToBarcode']
    
    # FILE HANDLES
    FASTQ = open(inp_dict['fastq_fp'], "r")
    TMPFNA = open(inp_dict['tmpFNA_fp'], "w")
    TRUNC = open(inp_dict['trunc_fp'], "w")
    

    # Current number of reads read
    nReads = 0
    # How many reads are candidates for mapping
    nTryToMap = 0
    # Current number of Long hits - meaning ...
    nLong = 0
    # Start and End of barcode within model:
    barcodeStart = inp_dict['model'].find("N")
    barcodeEnd = inp_dict['model'].rfind("N")
    barcodeLen = barcodeEnd - barcodeStart + 1
   
    # Dict of read name to 1 if mapped or pastEnd, 0 otherwise
    # only fills out when unmapped option
    mapnames = {}

    line_num = 0
    
    Find_cfg_d = {
            "flanking": inp_dict['flanking'],
            "wobbleAllowed": inp_dict['wobbleAllowed'],
            "minQuality": inp_dict['minQuality'],
            "fastq_fp": inp_dict['fastq_fp'],
            "debug": inp_dict['debug']
            }


    while (limit is None or (nReads < limit)):

        name = FASTQ.readline()
        line_num += 1


        if name == '':
            # If the file ends we break out of the loop
            break
        name = name.rstrip()
        if not name[0] == '@':
            raise Exception("Sequence name line does not start with @. File: "
            "{}, Line no. {}".format(inp_dict['fastq_fp'], line_num))


        seq = FASTQ.readline()
        line_num += 1
        seq = seq.rstrip()
        if not len(seq) > 0:
            raise Exception("Sequence line is empty. File {}, Line {}".format(
                inp_dict['fastq_fp'], line_num))
        if not re.match(r'^[A-Z]+$', seq):
            raise Exception("Sequence line contains invalid chars: " + seq \
                    + " \n File {}, Line {}".format(
                inp_dict['fastq_fp'], line_num) )


        line = FASTQ.readline()
        line_num += 1
        if not line[0] == '+':
            raise Exception("Third line does not start with +")

        quality = FASTQ.readline()
        line_num += 1
        quality = quality.rstrip()
        
        if not (quality != '' and len(quality) == len(seq)):
            raise Exception("Quality line is wrong length. "
                    " File {}, Line {}".format(
                inp_dict['fastq_fp'], line_num) )

        # Ignore second side of paired-end reads
        if re.match(r'^\S+ 2:', name):
            continue

        nReads += 1

        # Short sequences are unmappable
        if not len(seq) > len(inp_dict['model']) + inp_dict['minScore']:
            continue
        
        nLong += 1

        
        Find_cfg_d['line_num'] = line_num

        #obsStart is start of barcode within sequence
        # str, int
        barcode, obsStart = FindBarcode(seq, quality, inp_dict['model'],
                                        barcodeStart, barcodeEnd, 
                                        Find_cfg_d)
        

        if barcode is None or barcode == '':
            continue


        #logging.debug("start of barcode within sequence is: " + str(obsStart))

        # We create a shortname which subs out . to end of "name" of sequence
        shortname = re.sub(r' .*$', '', name)

        if shortname in nameToBarcode:
            raise Exception("Duplicate read name: {}\nFile {} line no. {}".format(
                            shortname, inp_dict['fastq_fp'], line_num -3))

        nameToBarcode[shortname] = barcode

        expModelStart = obsStart - barcodeStart
        transposonEnd = FindModelEnd(seq, inp_dict['model'], expModelStart,
                                    Find_cfg_d)
        if (transposonEnd is not None) and (len(seq) >= transposonEnd + inp_dict[
                'minScore']):
            # We write out to info files

            # The part of sequence in the Genome
            inGenome = seq[transposonEnd + 1:]
            TMPFNA.write(">{}\n{}\n".format(shortname, inGenome))
            nTryToMap += 1
            if inp_dict['unmapped_fp'] is not None:
                mapnames[shortname] = 0
            words = name.split(' ')
            words[0] += ":" + barcode
            TRUNC.write(" ".join(words) + "\n" + seq[transposonEnd:] \
                    + "\n+\n" + quality[transposonEnd:] + "\n")

    
    FASTQ.close()
    TMPFNA.close()
    TRUNC.close()
    logging.info("Read " + str(nReads) + " reads\n")
    
    if nTryToMap == 0:
        logging.warning("None of the reads are candidates for mapping."
                " None match the model and are long enough.\n")
        os.unlink(inp_dict['tmpFNA_fp'])
        return -1 

    ret_d = {
            "nReads": nReads,
            "nTryToMap": nTryToMap,
            "nLong": nLong,
            "mapnames": mapnames
            }

    return ret_d

    

def FindBarcode(seq, quality, model, expStart, expEnd, cfg_d):
    """
    Inputs (type and definition):

    seq (str) DNA sequence
    quality (str) quality of sequence
    model (str) model string
    expStart (int) Expected Start of barcode
    expEnd (int) Expected end of barcode
    cfg_d: (dict) (config_dict)
        wobbleAllowed (int) uncertainty in location of barcode or end of transposon,
            on either side of expectation
        minQuality (int) every nucleotide in a barcode must be at least this quality
        flanking (int) number of nucleotides on each side that must match
        debug (bool)
        fastq_fp: (str) The path to the FASTQ file
        line_num: (int) line number within FASTQ file of quality 

    Returns:
        barcode: (str) sequence of length ~20 which is the random barcode
        obsStart: (int) position of first nucleotide of barcode within seq
    """
    flanking = cfg_d['flanking']
    pre = model[expStart - flanking:expStart]
    post = model[expEnd + 1: expEnd + 1 + flanking]
    preLoc = FindSubstr(pre, seq, expStart - flanking, cfg_d['wobbleAllowed'])
    if preLoc is None:
        return [None, None]
    postLoc = FindSubstr(post, seq, expEnd + 1, cfg_d['wobbleAllowed'])
    if postLoc is None:
        return [None, None]

    # positions of 1st and last nucleotides of barcode
    start = preLoc + flanking
    end = postLoc - 1

    if not end - start == expEnd - expStart:
        if cfg_d['debug']:
            logging.info("Wrong barcode length: {} {} not {} {} in {}\n".format(
            start, end, expStart, expEnd, seq))
        return [None, None]

    barcode = seq[start:end + 1]

    if cfg_d['minQuality'] > 0:
        barqual = quality[start: end + 1]
        #quality score is encoded as 33
        scores = [ord(c) for c in barqual]
        for score in scores:
            if score < 0 or score > 100:
                raise Exception("Invalid score {} from barcode {}".format(
                score, barcode) + " quality {}".format(
                barqual) + " in file {} line {}".format(
                    cfg_d['fastq_fp'], cfg_d['line_num']))
            if score < cfg_d['minQuality']:
                if cfg_d['debug']:
                    logging.info("Low quality {} for barcode {} in {}\n".format(
                        score, barcode, seq) + "File {} line {}".format(
                            cfg_d['fastq_fp'], cfg_d['line_num']))
                return [None, None]
    
    return [barcode, start] # str, int


def FindSubstr(subseq, seq, expAt, wobble):
    """
    Here we search for model parts pre and post barcode within a sequence
    from reads.

    subseq (str) Subsequence within string we're looking for
    seq (str) Full sequence to search subseq within
    expAt (int) expected at location of subseq
    wobble (int) possible interval around expAt for which subseq can be found.
    """
    L = len(seq)
    K = len(subseq)
    for i in range(expAt - wobble, expAt + wobble + 1):
        if i >= 0 and i < L and seq[i:i+K] == subseq:
            return i

    return None



def FindModelEnd(seq, model, expOff, cfg_d):
    """
    Inputs: 

    seq: (str) String of Sequence from FASTQ
    model: (str) String of entire model
    expOff: (int) expected location model starts
    cfg_d:
        wobbleAllowed (int) uncertainty in location of barcode or end of transposon,
            on either side of expectation
        minQuality (int) every nucleotide in a barcode must be at least this quality
        flanking (int) number of nucleotides on each side that must match
        debug: (bool)
        fastq_fp: (str) The name of the FASTQ file
        line_num: (int) line number within FASTQ file of quality 

    Returns:
    None or transposonEnd (int) index in seq where transposon ends.

    """
    # expected End index
    expEnd = len(model) + expOff
    at = FindSubstr(model[-cfg_d['flanking']:],seq, expEnd - cfg_d['flanking'],
                    cfg_d['wobbleAllowed'])
    if at is None:
        if cfg_d['debug']:
            logging.info("No end sequence {} near position {} of \n{}\n".format(
                model[-cfg_d['flanking']:], expEnd - cfg_d['flanking'],
                seq) + " in file {} line no. {}".format(
                    cfg_d['fastq_fp'], cfg_d['line_num']))
        return None
    
    # last position of transposon
    return at + cfg_d['flanking'] - 1



        
            













def RunBLAT8(inp_dict):
    """
    Run this UCSC tool on a few cases:
    Case A:
        The queries file is full of end of sequences that had the  barcodes in
            them (tmpFNA sequences).
        The database is the pastEnd sequence
    Case B:
        The queries is the tmpFNA sequences - the barcodes reads
        The database is the genome FASTA file

    Inputs:
    
    inp_dict: (dict) must contain
        queriesFile (str) In case A & B: (tmpFNA_fp) 
        dbFile (str) case A: (endFNA_fp) caseB: (genome FASTA File)
        blat8_fp: (str) Path to output file
        blatcmd: (str) Path to BLAT command, i.e. Bin/blat
        minScore: (int) minimum score for mapping to genome or past-end
        minIdentity: (int) minimum %identity for mappinf to genome or past-end
        tileSize: (int) size of an alignment tile
        stepSize: (int) distance between the starting bases of alignment tiles
            (will overlap if stepSize<tileSize)
    

    """
    
    blat_args = [
            inp_dict['blatcmd'], 
            "-out=blast8", 
            "-t=dna",
            "-minScore=" + str(inp_dict['minScore']),
            "-minIdentity=" + str(inp_dict['minIdentity']),
            "-maxIntron=0",
            "-noTrimA",
            "-tileSize=" + str(inp_dict['tileSize']),
            "-stepSize=" + str(inp_dict['stepSize']),
            inp_dict['dbFile'],
            inp_dict['queriesFile'],
            inp_dict['blat8_fp']
            ]

    logging.info("Running blat8 with following command:\n{}".format(
        " ".join(blat_args)))
    blat8_result = subprocess.run(blat_args)
    logging.info("blat8 result: " + str(blat8_result))

    return None 


def HandleGenomeBLAT(rows, hitsPastEnd, HG_d, debug):
    """
    Given BLAT rows (as list of lists) and hitsPastEnd dict (as reference),
        output the mapping for the read, or not (if "trumped" by hit-past-end).
        Updates hitsPastEnd for the read to score=0 if 
            genomic hit trumps hit-past-end.


    Inputs:
    rows: list of lists, each sub list a blat8 split line
    hitsPastEnd: (dict)
    debug: (bool)
    HG_d (dict): Handle Genome Blat dict
        nPastEndTrumps: (int)  hit to past-end (close to) as good as hit to genome
        nPastEndIgnored: (int) weak hit to past-end ignored 
        nMapUnique: (int)
        nMapped: (int)
        nameToBarcode: (dict)
        delta: (int) minimum difference in score for considering a hit unique.
        OUTPUT: file handle for output

    """
    if len(rows) == 0:
        return None

    read = rows[0][0]

    if read == '' or read not in HG_d['nameToBarcode']:
        logging.info(rows)
        logging.info(HG_d['nMapped'])
        raise Exception("Handle Genome BLAT no Read")

    hits = []

    # indexes within besthits entries:
    SCAFFOLD, POSITION, STRAND, SCORE, QBEG, QEND = range(6)

    besthits = []

    for row in rows:
        read2, sbj, idn, l, mm, gaps, qB, qE, sB, sE, evl, score = row
        if not read2 == read:
            raise Exception("new read not equal: \n{}!=\n{}".format(read2, read))
        if len(besthits) == 0 or (
                float(score) >= float(besthits[0][SCORE]) - HG_d['delta']):
            # convert from 0-based to 1-based pos, and note that sB always < sE
            # so flip if stranded
            nv = "+" if int(sB) < int(sE) else "-"
            besthits.append([sbj, sB, nv, score, idn, qB, qE])
            if debug:
                logging.info("Hit for {}:{}:{} to {}:{}:{} score {}".format(
                    read, qB, qE, sbj, sB, sE, score) + " identity {}".format(
                        idn))

    if len(besthits) == 0:
        raise Exception("length of besthits is 0 within HandleGenomeBLAT...")

    scaffold, pos, strand, score, idn, qB, qE = besthits[0]

    # and output a mapping row (or none)
    if read in hitsPastEnd:
        if hitsPastEnd[read] >= float(score) - HG_d['delta']:
            HG_d['nPastEndTrumps'] += 1
            if debug:
                logging.info("Past end trumps for " + read)
            return None
        else:
            HG_d['nPastEndIgnored'] += 1
            if debug:
                logging.info("Ignoring hit past end of transposon for " + read)
            # This is so we don't print out a line for it later on...
            hitsPastEnd[read] = 0

    # else if no trumping, we write to output
    HG_d['OUTPUT'].write("\t".join([
        read, 
        HG_d['nameToBarcode'][read],
        scaffold,
        pos, 
        strand,
        "1" if len(besthits) == 1 else "0",
        qB,
        qE,
        score,
        idn + "\n"
        ]))

    if len(besthits) == 1:
        HG_d['nMapUnique'] += 1

    HG_d['nMapped'] += 1

    return None






def test():

    return None


def main():
    test()

    return None

if __name__ == "__main__":
    main()

