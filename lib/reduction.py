#! /usr/bin/env python

import glob
from os import path
import numpy as np
import time
try:
    import pyfits
except ImportError:
    from astropy.io import fits as pyfits


def make_master_dark(pathToDir):  # TODO What if no dark files found?
    """ Creates median of three newest darks in
    given directory"""
    outPath = path.join("workDir", "master_dark.fits")
    allDarks = glob.glob(path.join(pathToDir, "dark*.FIT"))
    if len(allDarks) < 3:
        # Not enough dark files were grabbed
        return None, 0, None
    lastDarks = sorted(allDarks, key=path.getctime)[-3:]
    darkNumber = path.split(lastDarks[0])[-1][4:-7]
    if (not path.exists(outPath)) or (path.getctime(outPath) < path.getctime(lastDarks[0])):
        # if there is no master dark at all or there are newer dark frames
        # we need to create a new one
        masterDarkData = np.median([safe_open_fits(fName)[0].data for fName in lastDarks], axis=0)
        masterDarkHDU = pyfits.PrimaryHDU(data=masterDarkData)
        masterDarkHDU.writeto(outPath, overwrite=True)
    else:
        # else hust use existing master dark
        masterDarkData = safe_open_fits(outPath)[0].data.copy()
    medianDarkValue = np.median(masterDarkData)
    stdDarkValue = np.std(masterDarkData)
    hotPixels = np.where(np.abs(masterDarkData-medianDarkValue) > stdDarkValue*10)
    return masterDarkData, darkNumber, hotPixels


def make_master_bias(pathToDir):  # the same quesion
    """ Creates median of all bias files in the given directory"""
    outPath = path.join("workDir", "master_bias.fits")
    if not path.exists(outPath):
        allBiases = glob.glob(path.join(pathToDir, "bias*.FIT"))
        if len(allBiases) != 5:
            return None
        masterBiasData = np.median([safe_open_fits(fName)[0].data for fName in allBiases], axis=0)
        masterBiasHDU = pyfits.PrimaryHDU(data=masterBiasData)
        masterBiasHDU.writeto(outPath, overwrite=True)
    else:
        masterBiasData = safe_open_fits(outPath)[0].data.copy()
    return masterBiasData


def safe_open_fits(pathToFile):
    """Occasionally we can try to read data from the file
        that is being written at the moment by CCDops.
        In this case we will get IOError and we have to
        wait for writing of the file to be finished"""
    for i in range(25):
        try:
            hdu = pyfits.open(pathToFile)
            return hdu
        except IOError:
            print("Got troubles while reading %s" % (pathToFile))
            time.sleep(0.1)
    print("Could not open file in 25 attempts")


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


def parse_object_file_name(fName):
    """ Function gets object name, add string (if exists) and
    filter name out of object file name"""
    print("parsing %s" % (fName))
    # File name is like objname+addstr+filter+number.FIT,
    # for example for wcom it may be 'wcom2b001.FIT'
    fNameNoExt = path.splitext(fName)[0]
    fNameNoNumbers = fNameNoExt[:-3]
    frameNumber = int(fNameNoExt[-3:])
    filtName = fNameNoNumbers[-1]
    fNameNoFilt = fNameNoNumbers[:-1]
    # Let's find what object is it
    for line in open(path.join("references", "list_of_objects.dat")):
        objName = line.strip()
        if fNameNoFilt.startswith(objName):
            addString = fNameNoFilt[len(objName):]
            return objName, filtName, addString, frameNumber
    # No sutable object found
    return None, None, None, None
