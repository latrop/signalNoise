#! /usr/bin/env python

from math import hypot, log10, atan2, degrees, radians
import subprocess
from shutil import move
import os
from os import path
from os import remove
import numpy as np

import pyfits


class SExCatalogue(object):
    def __init__(self, catFileName):
        self.objectList = []
        self.legend = []
        self.inds = [-1]
        for line in open(catFileName):
            sLine = line.strip()
            obj = {}
            if sLine.startswith("#"):
                # legend line
                self.legend.append(sLine.split()[2])
                self.inds.insert(-1, int(sLine.split()[1]) - 1)
                continue
            params = [float(p) for p in line.split()]
            for i in xrange(len(self.legend)):
                b = self.inds[i]
                e = self.inds[i+1]
                if e == b + 1:
                    obj[self.legend[i]] = params[b]
                else:
                    obj[self.legend[i]] = params[b:e]
            self.objectList.append(obj)
    
    def find_nearest(self, x, y):
        """ Returns nearest object to given coordinates"""
        nearest = min(self.objectList, key=lambda obj: hypot(x-obj["X_IMAGE"], y-obj["Y_IMAGE"]))
        dist = hypot(x-nearest["X_IMAGE"], y-nearest["Y_IMAGE"])
        if dist < 5.0:
            return nearest
        else:
            return None

    def get_median_value(self, parameter):
        return np.median([obj[parameter] for obj in self.objectList])


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
        fluxAuto = stObs["seParams"]["FLUX_AUTO"]
        fluxAper1 = stObs["seParams"]["FLUX_APER"][0]
        fluxAper2 = stObs["seParams"]["FLUX_APER"][1]
        fluxErr = abs(fluxAper2 - fluxAper1)
        stSn[st['name']] = fluxAper1/(fluxAper1+fluxErr)**0.5
        if not st["mag%s"%filtName.lower()] is None:
            magzpt = 2.5*log10(fluxAuto) + st["mag%s"%filtName.lower()]
            fluxzpt.append(10**(0.4*(30.0-magzpt)))

    # find object magnitude and s/n ratio
    if ref.objSEParams is None:
        return None, None, stSn
    xCen = ref.objSEParams["X_IMAGE"]
    yCen = ref.objSEParams["Y_IMAGE"]
    objPhotParams = cat.find_nearest(xCen, yCen)
    objFluxAuto = objPhotParams["FLUX_AUTO"]
    objFluxAper1 = objPhotParams["FLUX_APER"][0]
    objFluxAper2 = objPhotParams["FLUX_APER"][1]
    objFluxErr = abs(objFluxAper2 - objFluxAper1)
    objSn = objFluxAper1/(objFluxAper1+objFluxErr)**0.5
    if not st["mag%s"%filtName.lower()] is None:
        meanZpt = 30 - 2.5*log10(np.mean(fluxzpt))
        objMag = -2.5*log10(objFluxAuto)+meanZpt
    else:
        objMag = -99.0

    return objSn, objMag, stSn


def get_photometry_polar_mode(cat, ref):
    # we are not going to compute any magnitudes in polar mode
    # 1) find sn ratios for object pair
    fluxRatios = {}
    if not ref.objSEParams is None:
        xObjCen = ref.objSEParams["X_IMAGE"]
        yObjCen = ref.objSEParams["Y_IMAGE"]
        objPhotParams = cat.find_nearest(xObjCen, yObjCen)
        objFluxAuto = objPhotParams["FLUX_AUTO"]
        objFluxAper1 = objPhotParams["FLUX_APER"][0]
        objFluxAper2 = objPhotParams["FLUX_APER"][1]
        objFluxErr = abs(objFluxAper2-objFluxAper1)
        objSN = objFluxAper1 / (objFluxAper1+objFluxErr)**0.5
    else:
        objSN = None
        objFluxAuto = None
    # pair
    if not ref.objPairSEParams is None:
        xObjPairCen = ref.objPairSEParams["X_IMAGE"]
        yObjPairCen = ref.objPairSEParams["Y_IMAGE"]
        objPairPhotParams = cat.find_nearest(xObjPairCen, yObjPairCen)
        objPairFluxAuto = objPairPhotParams["FLUX_AUTO"]
        objPairFluxAper1 = objPairPhotParams["FLUX_APER"][0]
        objPairFluxAper2 = objPairPhotParams["FLUX_APER"][1]
        objPairFluxErr = abs(objPairFluxAper2 - objPairFluxAper1)
        objPairSN = objPairFluxAper1 / (objPairFluxAper1+objPairFluxErr)**0.5
    else:
        objPairSN = None
        objPairFluxAuto = None
    if (objFluxAuto is not None) and (objPairFluxAuto is not None):
        fluxRatios["obj"] = objFluxAuto / objPairFluxAuto
    else:
        fluxRatios["obj"] = None
    # 2) find sn ratios for standart pairs
    stSnDict = {}
    stFluxRatios = []
    for st, stPair in zip(ref.standartsObs, ref.standartPairsObs):
        if not st["seParams"] is None:
            stFluxAuto = st["seParams"]["FLUX_AUTO"]
            stFluxAper1 = st["seParams"]["FLUX_APER"][0]
            stFluxAper2 = st["seParams"]["FLUX_APER"][1]
            stFluxErr = abs(stFluxAper2 - stFluxAper1)
            stSN = stFluxAper1/(stFluxAper1+stFluxErr)**0.5
        else:
            stSN = None
            stFluxAuto = None
        if not stPair["seParams"] is None:
            stPairFluxAuto = stPair["seParams"]["FLUX_AUTO"]
            stPairFluxAper1 = stPair["seParams"]["FLUX_APER"][0]
            stPairFluxAper2 = stPair["seParams"]["FLUX_APER"][1]
            stPairFluxErr = abs(stPairFluxAper2 - stPairFluxAper1)
            stPairSN = stPairFluxAper1/(stPairFluxAper1+stPairFluxErr)**0.5
        else:
            stPairSN = None
            stPairFluxAuto = None
        stSnDict[st["name"]] = np.array([stSN, stPairSN])
        if (stFluxAuto is not None) and (stPairFluxAuto is not None):
            stFluxRatios.append(stFluxAuto/stPairFluxAuto)
    if stFluxRatios:
        fluxRatios["st"] = np.mean(stFluxRatios)
    else:
        fluxRatios["st"] = None
    return objSN, objPairSN, stSnDict, fluxRatios


def all_st_sn_decreased(stSnOld, stSnNew, polarMode):
    if polarMode:
        for key in stSnOld:
            if (key in stSnNew):
                a1 = stSnOld[key]
                sOld = sum(a1[np.where(a1 != np.array(None))])
                a2 = stSnNew[key]
                sNew = sum(a2[np.where(a2 != np.array(None))])
                if (sOld < sNew):
                    return False
            else:
                print key, stSnOld, stSnNew
    else:
        for key in stSnOld:
            if (key in stSnNew):
                sOld = stSnOld[key]
                sNew = stSnNew[key]
                if (sOld is not None) and (sNew is not None) and (sOld < sNew):
                    return False
    return True
