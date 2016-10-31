#! /usr/bin/env python

import time
import glob
import sys
import os
from os import path
sys.path.append(path.join(os.getcwd(), "lib"))
import numpy as np
from shutil import move
from collections import OrderedDict

import pylab

try:
    import pyfits
except ImportError:
    from astropy.io import fits as pyfits
    
if os.name == "nt":
    import winsound

from lib.GUIlib import *
from lib.reduction import *
from lib.se import *
from lib.filterPolCat import filterPolCat


def parse_object_file_name(fName):
    """ Function gets object name, add string (if exists) and 
    filter name out of object file name"""
    print "parsing %s" % (fName)
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


class MainApplication(Tk.Frame):
    def __init__(self, *args, **kwargs):
        self.objName = None
        self.currentObject = None
        self.currentFilter = None
        self.rawImages = []
        self.magData = []
        self.objSn = 0
        self.stSn = {}
        self.filterChecked = False
        self.desiredExposures = 0
        self.badObjects = []
        self.biasValue, self.darkValue = 0.0, 0.0
        self.photoLog = OrderedDict()
        # check if there is the working directory and create if nessesary
        if not path.exists("workDir"):
            os.mkdir("workDir")
        else:
            # cleanup working directory
            for f in glob.glob(path.join("workDir", "*.*")):
                os.remove(f)
        # cache flats
        self.flats = {}
        for filt in 'bvriXY':
            pathToFlat = path.join("tmp2", "flats", "ff%s.fts" % (filt))  # FIXME
            hdu = pyfits.open(pathToFlat)
            fData = np.flipud(hdu[0].data)
            self.flats[filt.lower()] = fData / np.mean(fData)
            hdu.close()
        # GUI stuff:
        self.root = Tk.Tk()
        self.root.title("SignalNoise")
        self.root.protocol('WM_DELETE_WINDOW', self.on_closing)
        self.showHotPixels = Tk.BooleanVar()
        self.showHotPixels.set(False)
        self.showHotPixels.trace("w", lambda a,b,c: self.update_plot())
        self.menubar = MenuBar(self)
        self.imagPanel = ImagPanel(self)
        self.rightPanel = RightPanel(self)
        self.root.mainloop()

    def on_closing(self):
        self.root.destroy()

    def setup(self, dirName):
        self.dirName = dirName
        self.cycle()

    def cycle(self):
        self.check_out_files()
        if not self.objName is None:
            self.rightPanel.update_object_info(self.objName, self.filtName, self.addString)
            self.run_computation()
        else:
            self.rightPanel.update_object_info("---", "---", "---")
            self.rightPanel.photometryString.set("")
        self.root.after(1500, self.cycle)        

    def check_out_files(self):
        allFiles = glob.glob(path.join(self.dirName, "*.FIT"))
        lightFiles = []
        # get rid of dark files and subdirectories
        for f in allFiles:
            fName = path.basename(f)
            if (not "dark" in fName) and (not "bias" in fName) and (path.isfile(f)):
                lightFiles.append(f)
        if len(lightFiles) == 0:
            return
        newestFile = max(lightFiles, key=path.getctime)
        fNameWithoutPath = path.basename(newestFile)
        # Let's find what object is it
        self.objName, self.filtName, self.addString, self.frameNumber = parse_object_file_name(fNameWithoutPath)

        # Determine if exposure is in polar mode
        if (self.filtName is not None) and (self.filtName.lower() in ("x", "y")):
            self.polarMode = True
        else:
            self.polarMode = False

    def clean_work_dir(self):
        for fName in glob.glob(path.join("workDir", "%s*_affineremap.fits" % self.currentObject)):
            os.remove(fName)
        for fName in glob.glob(path.join("workDir", "%s*.FIT" % self.currentObject)):
            os.remove(fName)

    def reset_new_filter(self):
        self.clean_work_dir()
        self.rawImages = []
        self.darkCleanImages = []
        self.magData = []
        self.currentFilter = self.filtName
        self.filterChecked = False
        self.objSn = 0
        self.stSn = {}
        if len(self.badObjects) > 0:
            self.rightPanel.update_bad_images_info([])
            self.badObjects = []

    def reset_new_object(self):
        self.reset_new_filter()
        self.currentObject = self.objName
        self.currentAddString = self.addString
        self.ref = Reference(self.objName)
        self.masterDarkData, darkNumber, self.hotPixels = make_master_dark(self.dirName)
        self.rightPanel.update_message("Dark", "Number %s" % darkNumber)
        self.masterBiasData = make_master_bias(self.dirName)
        objStr = "%s:%s"%(self.objName,self.addString)
        if objStr not in self.photoLog:
            # we only want to add a new object if there is no such
            # object in the dictionary. Otherwise it will erase log data
            # for the current onject if "reset" was called by user.
            self.photoLog[objStr] = {}

        # Clear working directory
        for f in glob.glob(path.join("workDir", "*")):
            if (not "dark" in path.basename(f)) and (not "bias" in path.basename(f)):
                os.remove(f)

    def rename_files(self, numberOfDesiredExposures):
        """
        Function renames files in such a way that there shall no be
        any gaps in the file numbers and numeration starts from the
        numberOfDesiredExposures value. i.e 002,005,006,008 become
        010,009,008,007 if numberOfDesiredExposure==10.
        """
        nexp = numberOfDesiredExposures
        tmpDir = path.join("workDir", "tmpRename")
        if not path.exists(tmpDir):
            os.mkdir(tmpDir)
        for fName in sorted(self.rawImages, reverse=True):
            fNameWithoutPath = path.basename(fName)
            objName, filtName, addString, frameNumber = parse_object_file_name(fNameWithoutPath)
            fout = "".join([objName, addString, filtName])
            fout += "%03i" % nexp
            nexp -= 1
            fout += ".FIT"
            fout = path.join(tmpDir, fout)
            move(fName, fout)
        for fName in glob.glob(path.join(tmpDir, "*.FIT")):
            fNameNoPath = path.basename(fName)
            move(fName, path.join(self.dirName, fNameNoPath))
        os.rmdir(tmpDir)
        return nexp


    def update_plot(self):
        self.imagPanel.plot_objects(self.ref, self.polarMode, self.hotPixels)
                           

    def run_computation(self):
        """ This is the main function. It is been
        called as soon as the object is selected."""
        if self.filtName != self.currentFilter:
            # new filter is being processed
            self.reset_new_filter()

        if (self.objName != self.currentObject) or (self.addString != self.currentAddString):
            # new object is being processed
            self.reset_new_object()

        # Observer can delete some raw files if they are bad,
        # so we need to check if there are files in rawImages list
        # that are not present on HDD anymore and delete them from
        # list and from workDir directory
        imageWasRemoved = False
        for f in self.rawImages:
            if not path.exists(f):
                imageWasRemoved = True
                self.rawImages.remove(f)
                fName = path.splitext(path.basename(f))[0]
                pathToFile = path.join("workDir", "%s_affineremap.fits" % fName)
                if path.exists(pathToFile):
                    os.remove(pathToFile)
                pathToFile = path.join("workDir", "%s.FIT" % fName)
                if path.exists(pathToFile):
                    os.remove(pathToFile)
                if pathToFile in self.darkCleanImages:
                    self.darkCleanImages.remove(pathToFile)
                if fName in self.badObjects:
                    self.badObjects.remove(fName)

        # Create the list of images to be processed
        newRawImages = []
        for img in sorted(glob.glob(path.join(self.dirName, "%s%s%s*" % (self.objName, self.addString, self.filtName)))):
            if img not in self.rawImages:
                newRawImages.append(img)
        if (len(newRawImages) == 0) and (not imageWasRemoved):
            print "no new images"
            return
        self.rawImages.extend(newRawImages)
        print(self.rawImages)

        # Check if filter name is equal to filter in file name
        if not self.filterChecked:
            hdu = safe_open_fits(self.rawImages[0])
            header = hdu[0].header
            headerFiltName = header["FILTER"].lower().strip()
            hdu.close()
            if headerFiltName != self.filtName.lower():
                # Filter missmatch alert
                res = self.rightPanel.update_message("Error", "Filter mismatch!")
                # res is True if the error message is new, so we enter
                # below if clause only one time every filter mismach event
                if (os.name == "nt") and res:
                    for i in xrange(5):
                        winsound.Beep(400, 500)
                self.rawImages = []  # Drop these bad files
                return
            else:
                self.rightPanel.update_message("Error", "")
                self.filterChecked = True

        # Check if the frame number is equal to desired number of frames (set by alarm window)
        if self.frameNumber == (self.desiredExposures-1):
            # one frame to go: one beep
            if (os.name == "nt"):
                winsound.Beep(700, 500)
        if self.frameNumber == self.desiredExposures:
            # final exposure: two beeps
            if (os.name == "nt"):
                winsound.Beep(700, 500)
                winsound.Beep(700, 500)
            self.desiredExposures *= -1

        # Subtract bias and dark files
        if newRawImages:
            newCleanImages, self.biasValue, self.darkValue = fix_for_bias_dark_and_flat(self.dirName, newRawImages,
                                                                                        self.masterDarkData,
                                                                                        self.masterBiasData,
                                                                                        self.flats[self.filtName.lower()])
            if not newCleanImages:
                return
            self.darkCleanImages.extend(newCleanImages)

        # Coadd images
        if not self.polarMode:
            self.numOfCoaddedImages = coadd_images(self.darkCleanImages, None)
        else:
            self.numOfCoaddedImages = coadd_images(self.darkCleanImages, self.filtName)
        self.rightPanel.update_message("Images summed", "%i" % self.numOfCoaddedImages)


        # Subtract background
        print "Cleaining background"
        backData = find_background(addString = " -BACKPHOTO_TYPE %s " % (self.menubar.backTypeVar.get()))
        summedFile = path.join("workDir", "summed.fits")
        self.imagPanel.plot_fits_file(summedFile)

        # Create a catalogue of objects on summed images
        catName = path.join("workDir", "field.cat")
        if self.polarMode:
            # in polar mode we have two catalogues: polar and filtred one.
            catNamePolar = path.join("workDir", "field_polar.cat")
            call_SE(summedFile, catNamePolar,
                    addString = "-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get()))
            filterPolCat(catNamePolar, catName, self.filtName)
            self.seCatPolar = SExCatalogue(catNamePolar)
            self.seCat = SExCatalogue(catName)

            # Find median FWHM of the image
            medianFWHM = self.seCat.get_median_value("FWHM_IMAGE")

            # Now we want to run SExtractor once again to get fluxes in
            # circular apertures of 1.55*FWHM and 1.55*sqrt(2)*FWHM radii
            aperRadius = 2*(1.55*medianFWHM+1)
            addString = "-PHOT_APERTURES %1.2f " % (2*aperRadius)
            addString += "-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get())
            call_SE(summedFile, catNamePolar, addString=addString)
            filterPolCat(catNamePolar, catName, self.filtName)
            self.seCatPolar = SExCatalogue(catNamePolar)
            self.seCat = SExCatalogue(catName)

            # Match objects from reference image on the observed one
            returnCode = self.ref.match_objects(summedFile, self.seCatPolar, self.filtName)
        else:
            call_SE(summedFile, catName,
                    addString="-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get()))
            self.seCat = SExCatalogue(catName)

            # Match objects from reference image on the observed one
            returnCode = self.ref.match_objects(summedFile, self.seCat)
            if not returnCode:
                # Find mean FWHM of standarts
                meanFWHM = self.ref.get_standatds_fwhm()
                print meanFWHM

                # Now we want to run SExtractor once again to get fluxes in
                # circular apertures of 1.55*FWHM and 1.55*sqrt(2)*FWHM radii
                aperRadius = (1.55*meanFWHM+1)
                addString = "-PHOT_APERTURES %1.2f " % (2*aperRadius)
                addString += "-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get())
                call_SE(summedFile, catName, addString=addString)
                self.seCat = SExCatalogue(catName)

                # Match objects once again. We already have a valid transform, so
                # we only need to match objects
                returnCode = self.ref.match_objects(summedFile, self.seCat, matchOnly=True)

        if returnCode:
            # matching failed. Maybe we need to coadd more_objects
            # to encrease S/N ratio.
            self.rightPanel.update_message("Error", "matching failed")
            # Remove object marks from image
            self.imagPanel.remove_objects_from_plot(upDateFigure=True)
            return
        else:
            self.rightPanel.update_message("Error", "")
        self.update_plot()

        # make aperture photometry of object and standarts
        snValues = {}
        if not self.polarMode:
            objSn, objMag, objMagSigma, stSn = get_photometry(self.seCat, self.ref, self.filtName, aperRadius,
                                                              self.biasValue, self.darkValue, backData)
            self.rightPanel.show_photometry_data(objSn, objMag, objMagSigma, stSn)
            self.magData.append((self.numOfCoaddedImages, objMag))
            # store magnitude value to a photomerty log
            self.photoLog["%s:%s"%(self.objName,self.addString)][self.filtName]=objMag

        elif self.polarMode:
            objSn, objPairSn, stSn, fluxRatios = get_photometry_polar_mode(self.seCatPolar, self.ref, aperRadius,
                                                                           self.biasValue, self.darkValue, backData)
            self.rightPanel.show_photometry_data_polar_mode(objSn, objPairSn, stSn, fluxRatios)

        # Check if object sn ratio decreased (for example due to a cloud)
        if ((objSn is not None) and (objSn < self.objSn)) and all_st_sn_decreased(self.stSn, stSn, self.polarMode):
            for f in newRawImages:
                objName = path.splitext(path.basename(f))[0]
                self.badObjects.append(objName)
        self.objSn = objSn
        self.stSn = stSn
        self.rightPanel.update_bad_images_info(self.badObjects)

        return


MainApplication()
