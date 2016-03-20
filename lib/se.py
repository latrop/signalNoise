#! /usr/bin/env python

from math import hypot, log10, atan2, degrees, radians
import subprocess
from shutil import move
import os
from os import path
from os import remove
import numpy as np

import pyfits


class SEObject(object):
    def __init__(self, seString):
        params = seString.split()
        self.number = int(params[0])
        self.xImage = float(params[1])
        self.yImage = float(params[2])
        self.fluxAuto = float(params[3])
        self.aImage = float(params[4])
        self.bImage = float(params[5])
        self.thetaImage = float(params[6])
        self.ellipticity = float(params[7])
        self.kronRadius = float(params[8])
        self.background = float(params[9])
        self.classStar = float(params[10])
        self.fwhmImage = float(params[11])
        self.fluxAper = float(params[12])
        self.fluxAutoErr = float(params[13])
        self.fluxAperErr = float(params[14])


class SExCatalogue(object):
    def __init__(self, catFileName):
        self.objectList = []
        for line in open(catFileName):
            sLine = line.strip()
            if sLine.startswith("#"):
                continue
            self.objectList.append(SEObject(sLine))
    
    def find_nearest(self, x, y):
        """ Returns nearest object to given coordinates"""
        nearest = min(self.objectList, key=lambda obj: hypot(x-obj.xImage, y-obj.yImage))
        dist = hypot(x-nearest.xImage, y-nearest.yImage)
        if dist < 3.0:
            return nearest
        else:
            return None


def call_SE(fitsFile, catName=None, addString=None):
    print "calling se"
    if os.name == "nt":  # Windows OS
        pathToSex = path.join('lib', 'sex.exe')
    else:
        pathToSex = 'sex'
    pathToFile = path.join("lib", "default.sex")
    callString = " ".join([pathToSex, fitsFile, '-c %s ' % pathToFile, "-VERBOSE_TYPE QUIET"])
    if not catName is None:
        callString += " -CATALOG_NAME %s " % catName
    if not addString is None:
        callString += addString
    print callString
    subprocess.call(callString, shell=True)
    print "done"


def clean_background():
    call_SE(path.join("workDir", "summed.fits"))
    move("background.fits", path.join("workDir", "background.fits"))
    origHDU = pyfits.open(path.join("workDir", "summed.fits"))
    origData = origHDU[0].data
    backHDU = pyfits.open(path.join("workDir", "background.fits"))
    backData = backHDU[0].data
    cleanData = origData - backData
    cleanHDU = pyfits.PrimaryHDU(data=cleanData)
    pathToFile = path.join("workDir", "back_clean.fits")
    if path.exists(pathToFile):
        remove(pathToFile)
    cleanHDU.writeto(pathToFile)    



def get_photometry(cat, ref, filtName):
    fluxzpt = []
    stSn = {}
    for st, stObs in zip(ref.standarts, ref.standartsObs):
        if stObs["seParams"] is None:
            stSn[st['name']] = None
            continue
        flux = stObs["seParams"].fluxAuto
        fluxErr = stObs["seParams"].fluxAutoErr
        stSn[st['name']] = abs(flux/fluxErr)
        if not st["mag%s"%filtName.lower()] is None:
            magzpt = 2.5*log10(flux) + st["mag%s"%filtName.lower()]
            fluxzpt.append(10**(0.4*(30.0-magzpt)))

    # find object magnitude and s/n ratio
    if ref.objSEParams is None:
        return None, None, stSn
    xCen = ref.objSEParams.xImage
    yCen = ref.objSEParams.yImage
    objPhotParams = cat.find_nearest(xCen, yCen)
    objFlux = objPhotParams.fluxAuto
    objFluxErr = objPhotParams.fluxAutoErr
    objSn = abs(objFlux/objFluxErr)
    if not st["mag%s"%filtName.lower()] is None:
        meanZpt = 30 - 2.5*log10(np.mean(fluxzpt))
        objMag = -2.5*log10(objFlux)+meanZpt
    else:
        objMag = -99.0

    return objSn, objMag, stSn


def get_photometry_polar_mode(cat, ref):
    # we are not going to compute any magnitudes in polar mode
    # 1) find sn ratios for object pair
    fluxRatios = {}
    if not ref.objSEParams is None:
        xObjCen = ref.objSEParams.xImage
        yObjCen = ref.objSEParams.yImage
        objPhotParams = cat.find_nearest(xObjCen, yObjCen)
        objFlux = objPhotParams.fluxAuto
        objFluxErr = objPhotParams.fluxAutoErr
        objSN = abs(objFlux / objFluxErr)
    else:
        objSN = None
        objFlux = None
    # pair
    if not ref.objPairSEParams is None:
        xObjPairCen = ref.objPairSEParams.xImage
        yObjPairCen = ref.objPairSEParams.yImage
        objPairPhotParams = cat.find_nearest(xObjPairCen, yObjPairCen)
        objPairFlux = objPairPhotParams.fluxAuto
        objPairFluxErr = objPairPhotParams.fluxAutoErr
        objPairSN = abs(objPairFlux / objPairFluxErr)
    else:
        objPairSN = None
        objPairFlux = None
    if (objFlux is not None) and (objPairFlux is not None):
        fluxRatios["obj"] = objFlux / objPairFlux
    else:
        fluxRatios["obj"] = None
    # 2) find sn ratios for standart pairs
    stSnDict = {}
    stFluxRatios = []
    for st, stPair in zip(ref.standartsObs, ref.standartPairsObs):
        if not st["seParams"] is None:
            stFlux = st["seParams"].fluxAuto
            stFluxErr = st["seParams"].fluxAutoErr
            stSN = abs(stFlux/stFluxErr)
        else:
            stSN = None
            stFlux = None
        if not stPair["seParams"] is None:
            stPairFlux = stPair["seParams"].fluxAuto
            stPairFluxErr = stPair["seParams"].fluxAutoErr
            stPairSN = abs(stPairFlux/stPairFluxErr)
        else:
            stPairSN = None
            stPairFlux = None
        stSnDict[st["name"]] = np.array([stSN, stPairSN])
        if (stFlux is not None) and (stPairFlux is not None):
            stFluxRatios.append(stFlux/stPairFlux)
    if stFluxRatios:
        fluxRatios["st"] = np.mean(stFluxRatios)
    else:
        fluxRatios["st"] = None
    return objSN, objPairSN, stSnDict, fluxRatios


def st_sn_increased(stSnOld, stSnNew, polarMode):
    if polarMode:
        for key in stSnOld:
            if (key in stSnNew):
                a1 = stSnOld[key]
                sOld = sum(a1[np.where(a1 != np.array(None))])
                a2 = stSnNew[key]
                sNew = sum(a2[np.where(a2 != np.array(None))])
                if (sOld > sNew):
                    return False
            else:
                print key, stSnOld, stSnNew
        return True
    else:
        for key in stSnOld:
            if (key in stSnNew):
                sOld = stSnOld[key]
                sNew = stSnNew[key]
                if (sOld is not None) and (sNew is not None) and (sOld > sNew):
                    return False
        return True
