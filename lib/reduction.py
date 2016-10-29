#! /usr/bin/env python

import glob
import os
from os import path
import numpy as np
import time
try:
    import pyfits
except ImportError:
    from astropy.io import fits as pyfits

    
def make_master_dark(pathToDir): # TODO What if no dark files found?
    """ Creates median of three newest darks in
    given directory"""
    outPath = path.join("workDir", "master_dark.fits")
    allDarks = glob.glob(path.join(pathToDir, "dark*.FIT"))
    lastDarks = sorted(allDarks, key=path.getctime)[-3:]
    if (not path.exists(outPath)) or (path.getctime(outPath)<path.getctime(lastDarks[0])):
        # if there is no master dark at all or there are newer dark frames
        # we need to create a new one
        masterDarkData = np.median([pyfits.open(fName)[0].data for fName in lastDarks], axis=0)
        masterDarkHDU =pyfits.PrimaryHDU(data = masterDarkData)
        if path.exists(outPath):
            os.remove(outPath)
        masterDarkHDU.writeto(outPath)
        darkNumber = path.split(lastDarks[0])[-1][4:-7]
    else:
        # else hust use existing master dark
        masterDarkData = pyfits.open(outPath)[0].data.copy()

    return masterDarkData, darkNumber


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


def safe_open_fits(pathToFile):
    """Occasionally we can try to read data from the file
        that is being written at the moment by CCDops.
        In this case we will get IOError and we have to
        wait for writing of the file to be finished"""
    for i in xrange(25):
        try:
            hdu = pyfits.open(pathToFile)
            return hdu
        except IOError:
            print "Got troubles while reading %s" % (pathToFile)
            time.sleep(0.1)
    print "Could not open file in 25 attempts"


def fix_for_bias_dark_and_flat(pathToDir, rawImages, biasData, darkData, flat):
    """ Function subtracts biases and darks from raw file """
    outFileList = []
    for pathToFile in rawImages:
        fName = path.basename(pathToFile)
        outName = path.join("workDir", fName)
        if path.exists(outName):
            continue
        hduRaw = safe_open_fits(pathToFile)
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

