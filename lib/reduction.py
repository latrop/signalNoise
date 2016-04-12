#! /usr/bin/env python

import glob
import os
from os import path
import numpy as np
import time
import pyfits

def make_master_dark(pathToDir): # TODO What if no dark files found?
    """ Creates median of three newest darks in
    given directory"""
    outPath = path.join("workDir", "master_dark.fits")
    if not path.exists(outPath):
        allDarks = glob.glob(path.join(pathToDir, "dark*.FIT"))
        lastDarks = sorted(allDarks, key=path.getctime)[-3:]
        masterDarkData = np.median([pyfits.open(fName)[0].data for fName in lastDarks], axis=0)
        masterDarkHDU =pyfits.PrimaryHDU(data = masterDarkData)
        if path.exists(outPath):
            os.remove(outPath)
        masterDarkHDU.writeto(outPath)
    else:
        masterDarkData = pyfits.open(outPath)[0].data.copy()
    return masterDarkData


def make_master_bias(pathToDir): # the same quesion
    """ Creates median of all bias files in the given directory"""
    outPath = path.join("workDir", "master_bias.fits")
    if not path.exists(outPath):
        allBiases = glob.glob(path.join(pathToDir, "bias*.FIT"))
        masterBiasData = np.median([pyfits.open(fName)[0].data for fName in allBiases], axis=0)
        masterBiasHDU =pyfits.PrimaryHDU(data = masterBiasData)
        if path.exists(outPath):
            os.remove(outPath)
        masterBiasHDU.writeto(outPath)
    else:
        masterBiasData = pyfits.open(outPath)[0].data.copy()
    return masterBiasData


def fix_for_bias_dark_and_flat(pathToDir, rawImages, flat):
    """ Function subtracts biases and darks from raw file """
    biasData = make_master_bias(pathToDir)
    darkData = make_master_dark(pathToDir)
    outFileList = []
    for pathToFile in rawImages:
        fName = path.basename(pathToFile)
        outName = path.join("workDir", fName)
        if path.exists(outName):
            continue
        for i in xrange(25):
            # Occasionally we can try to read data from the file
            # that is being written at the moment by CCDops.
            # In this case we will get IOError and we have to
            # wait for writing of the file to be finished
            try:
                hduRaw = pyfits.open(pathToFile)
                break
            except IOError:
                print "Got troubles while reading %s" % (pathToFile)
                time.sleep(0.1)
                continue
        else:
            # We didn't manage to read this file and we give up
            continue
        dataRaw = hduRaw[0].data.copy()
        headerRaw = hduRaw[0].header
        exptime = float(headerRaw['exptime'])
        darkDataReduced = darkData * (exptime/60.0)
        dataClr = (dataRaw/flat) - biasData - darkDataReduced
        outHDU = pyfits.PrimaryHDU(data=dataClr)
        outHDU.writeto(outName)
        hduRaw.close()
        outFileList.append(outName)
    if outFileList:
        return outFileList, np.median(biasData), np.median(darkDataReduced)
    else:
        return [], None, None

