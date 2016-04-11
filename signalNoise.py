#! /usr/bin/env python

import glob
import sys
import os
from os import path
sys.path.append(path.join(os.getcwd(), "lib"))
import numpy as np

import pylab

import pyfits

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
        self.currentObject = None
        self.currentFilter = None
        self.rawImages = []
        self.magData = []
        self.objSn = 0
        self.stSn = {}
        self.filterChecked = False
        self.desiredExposures = -1
        self.badObjects = []
        self.biasValue, self.darkValue = 0.0, 0.0
        # cache flats
        self.flats = {}
        for filt in 'bvriXY':
            pathToFlat = path.join("tmp2", "flats", "ff%s.fts" % (filt))  # FIXME
            hdu = pyfits.open(pathToFlat)
            fData = hdu[0].data
            self.flats[filt.lower()] = fData / np.mean(fData)
            hdu.close()
        # GUI stuff:
        self.root = Tk.Tk()
        self.root.title("SignalNoise")
#        self.root.geometry(("800x625"))
        self.root.protocol('WM_DELETE_WINDOW', self.on_closing)
        self.menubar = MenuBar(self)
        self.imagPanel = ImagPanel(self)
        # self.graphPanel = GraphPanel(self)
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
        newestFile = max(lightFiles, key=path.getctime)
        fNameWithoutPath = path.basename(newestFile)
        # Let's find what object is it
        self.objName, self.filtName, self.addString, self.frameNumber = parse_object_file_name(fNameWithoutPath)

        # Check if the frame number is equal to desired number of frames (set by alarm window)
        if self.frameNumber == (self.desiredExposures-1):
            # one frame to go: one beep
            if (os.name == "nt") and res:
                winsound.Beep(700, 500)
        if self.frameNumber == self.desiredExposures:
            # final exposure: two beeps
            if (os.name == "nt") and res:
                winsound.Beep(700, 500)
                winsound.Beep(700, 500)
            self.desiredExposures = -1

        # Determine if exposure is in polar mode
        if (self.filtName is not None) and (self.filtName.lower() in ("x", "y")):
            self.polarMode = True
        else:
            self.polarMode = False

    def reset_new_filter(self):
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
        self.ref = Reference(self.objName)

        # Clear working directory
        for f in glob.glob(path.join("workDir", "*")):
            if (not "dark" in path.basename(f)) and (not "bias" in path.basename(f)):
                os.remove(f)

    def run_computation(self):
        """ This is the main function. It is been
        called as soon as the object is selected."""
        if self.filtName != self.currentFilter:
            # new filter being processed
            self.reset_new_filter()

        if self.objName != self.currentObject:
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
            if not img in self.rawImages:
                newRawImages.append(img)
        if (len(newRawImages) == 0) and (not imageWasRemoved):
            print "no new images"
            return
            
        self.rawImages.extend(newRawImages)

        # Check if filter name is equal to filter in file name
        if not self.filterChecked:
            hdu = pyfits.open(self.rawImages[0])
            header = hdu[0].header
            headerFiltName = header["FILTER"].lower().strip()
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
            hdu.close()

        # Subtract bias and dark files
        if newRawImages:
            newCleanImages, self.biasValue, self.darkValue = fix_for_bias_dark_and_flat(self.dirName, newRawImages,
                                                                                        self.flats[self.filtName.lower()])
            self.darkCleanImages.extend(newCleanImages)
            self.rightPanel.update_message("Bias and dark", "Ok")

        # Coadd images
        if not self.polarMode:
            numOfCoaddedImages = coadd_images(self.darkCleanImages, None)
        else:
            numOfCoaddedImages = coadd_images(self.darkCleanImages, self.filtName)
        self.rightPanel.update_message("Images summed", "%i" % numOfCoaddedImages)


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

            # Find median FWHM of the image
            medianFWHM = self.seCat.get_median_value("FWHM_IMAGE")

            # Now we want to run SExtractor once again to get fluxes in
            # circular apertures of 1.55*FWHM and 1.55*sqrt(2)*FWHM radii
            aperRadius = (1.55*medianFWHM+1)
            addString = "-PHOT_APERTURES %1.2f " % (2*aperRadius)
            addString += "-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get())
            call_SE(summedFile, catName, addString=addString)
            self.seCat = SExCatalogue(catName)

            # Match objects from reference image on the observed one
            returnCode = self.ref.match_objects(summedFile, self.seCat)

        if returnCode:
            # matching failed. Maybe we need to coadd more_objects
            # to encrease S/N ratio.
            self.rightPanel.update_message("Error", "matching failed")
            # Remove object marks from image
            self.imagPanel.remove_objects_from_plot(upDateFigure=True)
            return
        else:
            self.rightPanel.update_message("Error", "")
        self.imagPanel.plot_objects(self.ref, self.polarMode)

        # make aperture photometry of object and standarts
        snValues = {}
        if not self.polarMode:
            objSn, objMag, objMagSigma, stSn = get_photometry(self.seCat, self.ref, self.filtName, aperRadius,
                                                              self.biasValue, self.darkValue, backData)
            self.rightPanel.show_photometry_data(objSn, objMag, objMagSigma, stSn)
            self.magData.append((numOfCoaddedImages, objMag))

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


[os.remove(f) for f in glob.glob(path.join("workDir", "*.*"))]
MainApplication()
