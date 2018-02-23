#! /usr/bin/env python

import glob
import sys
import os
from os import path
import numpy as np
from shutil import move
from collections import OrderedDict
import Tkinter as Tk

try:
    import pyfits
except ImportError:
    from astropy.io import fits as pyfits

if os.name == "nt":
    import winsound

sys.path.append(path.join(os.getcwd(), "lib"))
from lib import GUIlib  # noqa
from lib import alignment  # noqa
from lib import reduction # noqa
from lib import se  # noqa
from lib.filterPolCat import filterPolCat  # noqa


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
        self.object_selected_manually = False
        self.photoLog = OrderedDict()
        # check if there is the working directory and create if nessesary
        if not path.exists("workDir"):
            os.mkdir("workDir")
        else:
            # cleanup working directory
            for f in glob.glob(path.join("workDir", "*.*")):
                try:
                    # File can be protected from removing (for example if it is
                    # opened by another application)
                    os.remove(f)
                except OSError:
                    pass
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
        self.showHotPixels.trace("w", lambda a, b, c: self.update_plot())
        self.menubar = GUIlib.MenuBar(self)
        self.imagPanel = GUIlib.ImagPanel(self)
        self.rightPanel = GUIlib.RightPanel(self)
        self.root.mainloop()

    def on_closing(self):
        self.root.destroy()

    def setup(self, dirName):
        self.dirName = dirName
        self.cycle()

    def cycle(self):
        self.check_out_files()
        if self.objName is not None:
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
            if ("dark" not in fName) and ("bias" not in fName) and (path.isfile(f)):
                lightFiles.append(f)
        if len(lightFiles) == 0:
            return
        newestFile = max(lightFiles, key=path.getctime)
        fNameWithoutPath = path.basename(newestFile)
        # Let's find what object is it
        if not self.object_selected_manually:
            (self.objName, self.filtName,
             self.addString, self.frameNumber) = reduction.parse_object_file_name(fNameWithoutPath)

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
        self.ref = alignment.Reference(self.objName)
        self.masterDarkData, darkNumber, self.hotPixels = reduction.make_master_dark(self.dirName)
        self.rightPanel.update_message("Dark", "Number %s" % darkNumber)
        self.masterBiasData = reduction.make_master_bias(self.dirName)
        objStr = "%s:%s" % (self.objName, self.addString)
        if objStr not in self.photoLog:
            # we only want to add a new object if there is no such
            # object in the dictionary. Otherwise it will erase log data
            # for the current onject if "reset" was called by user.
            self.photoLog[objStr] = {}

        # Clear working directory
        for f in glob.glob(path.join("workDir", "*")):
            if ("dark" not in path.basename(f)) and ("bias" not in path.basename(f)):
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
            objName, filtName, addString, frameNumber = reduction.parse_object_file_name(fNameWithoutPath)
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
        searchTemplate = path.join(self.dirName, "%s%s%s*" % (self.objName, self.addString, self.filtName))
        for img in sorted(glob.glob(searchTemplate)):
            if img not in self.rawImages:
                newRawImages.append(img)
        if (len(newRawImages) == 0) and (not imageWasRemoved):
            print("no new images")
            return
        self.rawImages.extend(newRawImages)
        print(self.rawImages)

        # Check if filter name is equal to filter in file name
        if not self.filterChecked:
            hdu = reduction.safe_open_fits(self.rawImages[0])
            header = hdu[0].header
            headerFiltName = header["FILTER"].lower().strip()
            hdu.close()
            if headerFiltName != self.filtName.lower():
                # Filter missmatch alert
                res = self.rightPanel.update_message("Error", "Filter mismatch!")
                # res is True if the error message is new, so we enter
                # below if clause only one time every filter mismach event
                if (os.name == "nt") and res:
                    for i in range(5):
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
            flat = self.flats[self.filtName.lower()]
            newCleanImages, biasValue, darkValue = reduction.fix_for_bias_dark_and_flat(self.dirName, newRawImages,
                                                                                        self.masterDarkData,
                                                                                        self.masterBiasData, flat)
            if not newCleanImages:
                return
            self.darkCleanImages.extend(newCleanImages)

        # Coadd images
        if not self.polarMode:
            self.numOfCoaddedImages = alignment.coadd_images(self.darkCleanImages, None)
        else:
            self.numOfCoaddedImages = alignment.coadd_images(self.darkCleanImages, self.filtName)
        self.rightPanel.update_message("Images summed", "%i" % self.numOfCoaddedImages)

        # Subtract background
        print("Cleaining background")
        backData = se.find_background(addString=" -BACKPHOTO_TYPE %s " % (self.menubar.backTypeVar.get()))
        summedFile = path.join("workDir", "summed.fits")
        self.imagPanel.plot_fits_file(summedFile)

        # Create a catalogue of objects on summed images
        catName = path.join("workDir", "field.cat")
        if self.polarMode:
            # in polar mode we have two catalogues: polar and filtred one.
            catNamePolar = path.join("workDir", "field_polar.cat")
            # if there is no given aperture size, we need to run the SExtractor twice:
            # one time to find FWHM and compute aperture, and one time to compute acutal
            # fluxes with this aperture
            if self.ref.apertureSize is None:
                se.call_SE(summedFile, catNamePolar,
                           addString="-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get()))
                filterPolCat(catNamePolar, catName, self.filtName)
                self.seCatPolar = se.SExCatalogue(catNamePolar)
                self.seCat = se.SExCatalogue(catName)

                # Find median FWHM of the image
                medianFWHM = self.seCat.get_median_value("FWHM_IMAGE")

                # Now we want to run SExtractor once again to get fluxes in
                # circular apertures of 1.55*FWHM
                aperRadius = 1.55*medianFWHM+1
            else:
                aperRadius = self.ref.apertureSize
            addString = "-PHOT_APERTURES %1.2f " % (2*aperRadius)
            addString += "-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get())
            se.call_SE(summedFile, catNamePolar, addString=addString)
            filterPolCat(catNamePolar, catName, self.filtName)
            self.seCatPolar = se.SExCatalogue(catNamePolar)
            self.seCat = se.SExCatalogue(catName)

            # Match objects from reference image on the observed one
            returnCode = self.ref.match_objects(summedFile, self.seCatPolar, self.filtName)
        else:
            se.call_SE(summedFile, catName,
                       addString="-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get()))
            self.seCat = se.SExCatalogue(catName)

            # Match objects from reference image on the observed one
            returnCode = self.ref.match_objects(summedFile, self.seCat)
            if not returnCode:
                # Find mean FWHM of standarts
                meanFWHM = self.ref.get_standatds_fwhm()
                print(meanFWHM)

                if self.ref.apertureSize is None:
                    # if the aperture size is not given, we want to compute one
                    print("using 1.55*FWHM+1 rule to compute an aperture")
                    aperRadius = 1.55*meanFWHM+1
                else:
                    print("++++++++ using given aperture")
                    aperRadius = self.ref.apertureSize
                # Now we want to run SExtractor once again to get fluxes in circular apertures
                addString = "-PHOT_APERTURES %1.2f " % (2*aperRadius)
                addString += "-BACKPHOTO_TYPE %s" % (self.menubar.backTypeVar.get())
                se.call_SE(summedFile, catName, addString=addString)
                self.seCat = se.SExCatalogue(catName)

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
        if not self.polarMode:
            objSn, objMag, objMagSigma, stSn = se.get_photometry(self.seCat, self.ref, self.filtName, aperRadius,
                                                                 biasValue, darkValue, backData)
            self.rightPanel.show_photometry_data(objSn, objMag, objMagSigma, stSn)
            self.magData.append((self.numOfCoaddedImages, objMag))
            # store magnitude value to a photomerty log
            self.photoLog["%s:%s" % (self.objName, self.addString)][self.filtName] = objMag

        elif self.polarMode:
            objSn, objPairSn, stSn, fluxRatios = se.get_photometry_polar_mode(self.seCatPolar, self.ref, aperRadius,
                                                                              biasValue, darkValue, backData)
            self.rightPanel.show_photometry_data_polar_mode(objSn, objPairSn, stSn, fluxRatios)

        # Check if object sn ratio decreased (for example due to a cloud)
        if ((objSn is not None) and (objSn < self.objSn)) and se.all_st_sn_decreased(self.stSn, stSn, self.polarMode):
            for f in newRawImages:
                objName = path.splitext(path.basename(f))[0]
                self.badObjects.append(objName)
        self.objSn = objSn
        self.stSn = stSn
        self.rightPanel.update_bad_images_info(self.badObjects)

        return


MainApplication()
