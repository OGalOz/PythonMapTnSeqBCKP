#!python3
# this file is a translation of DesignRandomPool.pl into python
# Given the output of MapTnSeq.pl, possibly from multiple runs,
# choose the reliable and unambiguous tags
# Makes a tab-delimited table (Pool File)

# Inputs are one or more MapTnSeq output tables that contain rows:
# read,barcode,scaffold,pos,strand,uniq,qBeg,qEnd

import os
import sys
import logging
import json
import copy
import subprocess



def RunDesignRandomPool(inp_d, DEBUGPRINT):
    """
    We run through the program, dividing each part into a separate function

    INPUTS:
    inp_d: (dict) 
        Required inputs listed in function "ParseInputs" below

    DEBUGPRINT: (bool) Normally False if not checking changes in variables.
        (If you want to print out JSON files in between each part of the function
            for debugging purposes)
    """

    parsed_vars = ParseInputs(inp_d)

    new_vars = InitNewVars1()
    parsed_vars.update(new_vars)

    if DEBUGPRINT:
        with open(os.path.join(parsed_vars['tmp_dir'], "Vars01.json")) as g:
            g.write(json.dumps(parsed_vars, indent=4))

    ProcessInputMapTnSeqTables(parsed_vars)


    if DEBUGPRINT:
        with open(os.path.join(parsed_vars['tmp_dir'], "Vars02.json")) as g:
            g.write(json.dumps(parsed_vars, indent=4))

    POOLFH = InitPoolFileHandle(parsed_vars["output_fp"])
    parsed_vars["POOL_FH"] = POOLFH

    new_vars = CountBarCodesPrintPool(parsed_vars)
    parsed_vars.update(new_vars)


    if DEBUGPRINT:
        with open(os.path.join(parsed_vars['tmp_dir'], "Vars03.json")) as g:
            g.write(json.dumps(parsed_vars, indent=4))

    GetVariantsPrintPool(parsed_vars)


    if DEBUGPRINT:
        with open(os.path.join(parsed_vars['tmp_dir'], "Vars04.json")) as g:
            g.write(json.dumps(parsed_vars, indent=4))


    getChao2Estimates(parsed_vars['totcodes'], 
            parsed_vars['f1'], 
            parsed_vars['f2'])
   
    RunPoolStatsR(parsed_vars)

    logging.info("Finished running Design Random Pool")



def ParseInputs(input_dict):

    """
    input_dict: (dict)
        output_fp: (str) Path to pool file
        genes_table_fp: (str) Path to genes_table
        R_fp: (str) Path to PoolStats.R
        tmp_dir: (str) Path to temporary directory
        minN: (int)
        minFrac: (float)
        minRatio: (float)
        maxQBeg: (float)
        map_tnseq_table_fps: (list of str)
    """
    
    # Check strings:
    for x in ["output_fp", "genes_table_fp", "R_fp", "tmp_dir"]:
        if x not in input_dict or not isinstance(input_dict[x],str):
            raise Exception(x + " must be an argument and must be string")

    # Check floats
    for x in ["minFrac", "minRatio", "maxQBeg"]:
        if x not in input_dict or not isinstance(input_dict[x],float):
            raise Exception(x + " must be an argument and must be float")
    # Check input files
    if "map_tnseq_table_fps" in input_dict:
        for fp in input_dict["map_tnseq_table_fps"]:
            if not os.path.exists(fp):
                raise Exception("Internal MapTnSeq file {} not found".format(
                    fp))
    else:
        raise Exception("map_tnseq_table_fps not found as input to Design Random Pool")

    # minN
    if "minN" not in input_dict or not isinstance(input_dict["minN"], int):
        raise Exception("minN input to Design RandomPool incorrect")

    logging.info("All input parameters to Design Random Pool passed")


def RunPoolStatsR(inp_d):
    """
    inp_d: (dict) contains
        R_fp (str)
        pool_fp (str)
        genes_table_fp (str)
        nMapped (int)
        R_op: (str) Path to R log?
    """

    R_executable = "Rscript"

    RCmds = [R_executable, 
            inp_d["R_fp"], 
            inp_d["pool_fp"],
            inp_d["genes_table_fp"],
            inp_d["nMapped"]]


    logging.info("Running R PoolStats")

    with open(inp_d["R_op"]) as f:
        subprocess.call(RCmds, stdout=f)

    if os.path.exists(inp_d["R_op"]):
        logging.info("Succesfully ran R, log found at " + inp_d["R_op"])







def getChao2Estimates(totcodes, f1, f2):
    """
    All inputs are int
    """
    chao = int(totcodes + f1**2/(2*f2 + 1))
    logging.info("Chao2 estimate of #barcodes present (may be inflated for " \
            + "sequencing error): " + str(chao))
    



def GetVariantsPrintPool(inp_dict):
    """
    inp_dict: (dict) Contains
        barcodeAt: (dict)
        POOL_FH: file handle to output pool file
    """
    
    nOut = 0
    nMasked = 0
    nMaskedReads = 0
    barcodeAt = inp_dict['barcodeAt']

    # Pool file will have 1 per non "Masked" keys in barcodeAt
    for k in barcodeAt.keys():
        barcode, row = k, barcodeAt[k]
        
        nTot, nMax, maxAt, nNext, nextAt = row

        # This returns a list of strings with each nucleotide
        # subbed with the other three nucleotides (e.g. A -> C,G,T)
        # If there are M barcodes, it return 3*M barcodes.
        variants = Variants(barcode)

        mask = 0

        for variant in variants:
            if (variant in barcodeAt \
                    and barcodeAt[variant][0] > nTot \
                    and barcodeAt[variant][2] == maxAt):
                nMasked += 1
                nMaskedReads += nMax
                mask = 1
                continue

        if mask != 0:
            continue
        
        # maxAt is a key "A:B:C"
        atSplit = maxAt.split(":")
       
        # nextAt could be "::" or "A:B:C"
        nextAtSplit = nextAt.split(":")

        if barcode in inp_dict['pastEnd']:
            nPastEnd = inp_dict['pastEnd'][barcode]
        else:
            nPastEnd = 0

        inp_dict['POOL_FH'].write("\t".join([
                        barcode,
                        ReverseComplement(barcode),
                        str(nTot),
                        str(nMax),
                        "".join(atSplit),
                        str(nNext),
                        "".join(nextAtSplit),
                        str(nPastEnd) + "\n"
                        ))
        nOut += 1

    
    inp_dict['POOL_FH'].close()

    logging.info("Masked {} off-by-1 barcodes ({} reads) leaving {}".format(
                 nMasked, nMaskedReads, nOut) + " barcodes.")
    logging.info("Reads for those barcodes: {} of {} ({})".format(
                nReadsForUsable, nMapped, 100*nReadsForUsable/(nMapped + 10**-6)))

        

            







def CountBarCodesPrintPool(inp_dict):
    """
    Inputs:
        inp_dict contains:
            barPosCount: (dict) A dict which maps to dicts with each subdict 
                being a key (a:b:c) which maps to a list with 
                [nReads, nGoodReads]
            pastEnd: (dict)
            minN: (int)
            POOL_FH: File Handle for pool file
            minFrac: (float)
            minRatio: (float)
    """
    
    # For barcodes, f1 is number seen once, f2 is number seen twice,
    # totcodes is number of barcodes seen.
    f1 = 0
    f2 = 0
    totcodes = 0

    # Classification of bcs
    nInCategory = {"Usable": 0, 
                   "PastEnd": 0, 
                   "NoCons":0, 
                   "FewGood": 0, 
                   "LoRatio": 0
                   }
    SCAFFOLD, POS, STRAND, UNIQ, QBEG, QEND, SCORE, IDENTITY = range(8)

    # Barcodes seen more than once
    nMulti = 0
    
    # barcode to list of nTot, nMax, at, nNext, nextAt
    barcodeAt= {}

    nReadsForUsable = 0 #?

    for k in inp_dict['barPosCount'].keys():
        barcode, key_d = k, inp_dict['barPosCount'][k]

        if barcode in inp_dict['pastEnd']:
            nPastEnd = inp_dict['pastEnd'][bc]
        else:
            nPastEnd = 0
        
        nTot = nPastEnd

        for j in key_d.keys():
            key, value = j, key_d[j]
            #nTot increases by the total reads
            nTot += value[0]

        if nTot == 1:
            f1 += 1
        elif nTot == 2:
            f2 += 1

        totcodes += 1

        if not nTot >= inp_dict['minN']:
            continue

        nMulti += 1

        # "at" will be a sorted list of keys based on their 'nReads' (0th index)
        # list is sorted by decreasing order (highest first)
        preAt = sorted(key_d.items(), key=lambda j: j[1][0], reverse=True)
        at = [v[0] for v in z] 
        if len(at) == 0:
            nMax = 0
        else:
            # nMax takes the maximum nReads value (Why?)
            nMax = key_d[at[0]][0]

        if float(nPastEnd) >= float(nTot)/2 or nPastEnd >= nMax:
            nInCategory["PastEnd"] += 1
            n2 = nMax
            inp_dict['POOL_FH'].write("\t".join(
                barcode,
                ReverseComplement(barcode),
                str(nTot),
                str(nPastEnd),
                "PastEnd",
                "",
                "",
                str(n2),
                "",
                "",
                "",
                str(nPastEnd) + "\n"
                ))
            continue

        if not (nMax >= minN and \
                float(float(nMax)/float(nTot)) >= inp_dict['minFrac']):
            nInCategory["NoCons"] += 1
            continue

        maxAt = at[0]

        # Checking unique & qbeg=1 -- but the latter part may be redundant with 
        # maxQBeg
        nGood = key_d[maxAt][1]

        if not (nGood >= minN and \
                float(float(nGood)/float(nTot)) >= inp_dict['minFrac']:
            nInCategory["FewGood"] += 1
            continue

        nextAt = "::"
        nNext = 0
        if len(at) > 1:
            nextAt = at[1]
            nNext = key_d[nextAt][0]

        if not (nMax >= inp_dict['minRatio'] * nNext):
            nInCategory["LoRatio"] += 1
            continue

        barcodeAt[barcode] = [nTot, nMax, maxAt, nNext, nextAt]
        nInCategory["Usable"] += 1
        nReadsForUsable += nTot



    logging.info("{} barcodes seen {} or more times, ".format(
                nMulti, inp_dict['minN']) \
            + " map {} (minFrac {} minRatio {})".format(
                nInCategory["Usable"], 
                inp_dict['minFrac'],
                inp_dict['minRatio']))

    for category in sorted(nInCategory.keys()):
        if nInCategory[category] > 0:
            logging.info("{}\t{}\t{}".format(
                category,
                nInCategory[category],
                nInCategory[category]/ nMulti
                ))


    ret_d = {
        
        "barcodeAt": barcodeAt,
        "f1": f1,
        "f2": f2,
        "totcodes": totcodes

        }


    return ret_d
















def InitPoolFileHandle(poolfile_path):
    # poolfile_path is a string, path to output

    logging.info("Starting to write to Pool File at " + poolfile_path)
    POOL_FH = open(poolfile_path, "w")

    POOL_FH.write("\t".join([
        "barcode",
        "rcbarcode",
        "nTot",
        "n",
        "scaffold",
        "strand",
        "pos",
        "n2",
        "scaffold2",
        "strand2",
        "pos2", 
        "nPastEnd\n",
        ]))

    return POOL_FH


def InitNewVars1():
    
    nMapped = 0 # reads considered
    nSkipQBeg = 0 # reads ignored because qBeg > maxQBeg
    barPosCount = {} # {barcode => {scaffold:strand:position => 
    #                                [nReads, nGoodReads]}} 
    #                   where a "good" read has uniq=1 and qBeg = 1

    pastEnd = {} # number of reads for a barcode mapped past the end 
    # (pastEnd reads are not included in barPosCount)

    new_vars = {
            "nMapped": nMapped,
            "nSkipQBeg": nSkipQBeg,
            "barPosCount": barPosCount,
            "pastEnd": pastEnd
            }

    return new_vars
    



def ProcessInputMapTnSeqTables(inp_dict):
    """
    In this function we run through all the input Map TnSeq tables
    We count the number of mapped, pastEnd, and normal.

    inp_dict: (dict)

        map_tnseq_table_fps: (list of strings (filepaths))
        pastEnd_d: (dict) A dictionary mapping barcode of pastEnds to counts
        maxQBeg: (int) 
        barPosCount: (dict) A dict which maps to dicts with each subdict 
            being a key (a:b:c) which maps to a list with [nReads, nGoodReads]
        nMapped: (int)
        nSkipQBeg: (int)
    """
    logging.info("Reading Mapping files:\n" \
            + "\n".join(inp_dict['map_tnseq_table_fps']))
    
    for MTS_fp in inp_dict['map_tnseq_table_fps']:
        MTS_FH = open(MTS_fp,"r")

        c_line = MTS_FH.readline()
        while c_line != '':
            c_line = c_line.rstrip()
            # read, barcode, scaffold, position, strand, unique, qBeg, qEnd, 
            #   score, identity
            read, bc, scf, pos, strnd, unq, qB, qE, scr, idn = c_line.split('\t')
            if scf == "pastEnd":
                if bc in inp_dict['pastEnd_d']:
                    inp_dict['pastEnd_d'][bc] += 1
                else:
                    inp_dict['pastEnd_d'][bc] = 1
            elif inp_dict['maxQBeg'] >= 1 and qB <= inp_dict['maxQBeg']:
                key = ":".join(scf, strnd, pos)
                UpdateBarPosCount(inp_dict['barPosCount'], bc, key, unq, qB) 
                inp_dict['nMapped'] += 1
            else:
                inp_dict['nSkipQBeg'] += 1

            c_line = MTS_FH.readline()

        MTS_FH.close()

    logging.info("Read {} mapped reads for {} distinct barcodes.".format(
                inp_dict['nMapped'],
                len(inp_dict['barPosCount'].keys())))
    logging.info("Skipped {} reads with qBeg > {}".format(
                inp_dict['nSkipQBeg'],
                inp_dict['maxQBeg']))




def UpdateBarPosCount(barPosCount, barcode, key, unique, qBeg):
    """
    barPosCount: (dict)
    barcode: (str)
    key: (str)
    unique: (str) string of int (0 or 1)
    qBeg: (str)
    """
    if barcode in barPosCount:
        if key in barPosCount[barcode]:
            barPosCount[barcode][key][0] += 1
            if unique == "1" and qBeg == "1":
                barPosCount[barcode[key][1] += 1
        else:
            if unique == "1" and qBeg == "1":
                barPosCount[barcode][key] = [1,1]
            else:
                barPosCount[barcode][key] = [1,0]
    else:
        if unique == "1" and qBeg == "1":
            barPosCount[barcode][key] = [1,1]
        else:
            barPosCount[barcode][key] = [1,0]




def ReverseComplement(barcode):
    """
    barcode: (str) just "A""C""T""G"
    We return reverse complement: ACCAGT -> ACTGGT
    """
    revc_bc = ""
    barcode = barcode.upper()
    transl_d = {"A":"T","C":"G","T":"A","G":"C"}
    for char in barcode:
        if char not in ["A","C","T","G"]:
            raise Exception("char {} not in ACTG as expected".format(char))
        else:
            revc_bc += (transl_d[char])

    return revc_bc[::-1]
            


def Variants(barcode):
    """
    barcode: (str) a string of a DNA sequence
    """
    out = []
    baseseq = barcode.upper()

    for i in range(len(baseseq)):
        pre = baseseq[0:i]
        char = baseseq[i]
        post = baseseq[i+1:]
        if not char in ["A","C","G","T"]:
            continue
        for c in ["A","C","G","T"]:
            if not char == c:
                out.append(pre + c + post)

    return out





def test():

    return None


def main():
    test()

    return None

if __name__ == "__main__":
    main()
