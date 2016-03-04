#! /usr/bin/env python

import glob
import os
from os import path
import numpy as np
import pyfits
import time

def make_master_dark(pathToDir): # TODO What if no dark files found?
    """ Creates median of three newest darks in
    given directory"""
    allDarks = glob.glob(path.join(pathToDir, "dark*.FIT"))
    lastDarks = sorted(allDarks, key=path.getctime)[-3:]
    masterDarkData = np.median([pyfits.open(fName)[0].data for fName in lastDarks], axis=0)
    masterDarkHDU =pyfits.PrimaryHDU(data = masterDarkData)
    outPath = path.join("workDir", "master_dark.fits")
    if path.exists(outPath):
        os.remove(outPath)
    masterDarkHDU.writeto(outPath)
    return masterDarkData


def make_master_bias(pathToDir): # the same quesion
    """ Creates median of all bias files in the given directory"""
    allBiases = glob.glob(path.join(pathToDir, "bias*.FIT"))
    masterBiasData = np.median([pyfits.open(fName)[0].data for fName in allBiases], axis=0)
    masterBiasHDU =pyfits.PrimaryHDU(data = masterBiasData)
    outPath = path.join("workDir", "master_bias.fits")
    if path.exists(outPath):
        os.remove(outPath)
    masterBiasHDU.writeto(outPath)
    return masterBiasData


def fix_for_bias_dark_and_flat(pathToDir, rawImages, flat):
    """ Function subtracts biases and darks from raw file """
    biasData = make_master_bias(pathToDir)
    darkData = make_master_dark(pathToDir)
    outFileList = []
    for pathToFile in rawImages:
        fName = path.basename(pathToFile)
        outName = path.join("workDir", fName)
        outFileList.append(outName)
        if path.exists(outName):
            continue
        hduRaw = pyfits.open(pathToFile)
        dataRaw = hduRaw[0].data
        headerRaw = hduRaw[0].header
        exptime = float(headerRaw['exptime'])
        dataClr = (dataRaw/flat) - biasData - darkData * (exptime/60.0)
        outHDU = pyfits.PrimaryHDU(data=dataClr)
        outHDU.writeto(outName)
    return outFileList

