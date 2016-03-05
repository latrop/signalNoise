#! /usr/bin/env python

import glob
import sys
import os
from os import path
sys.path.append(path.join(os.getcwd(), "lib"))
import numpy as np

import pylab

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
    filtName = fNameNoNumbers[-1]
    fNameNoFilt = fNameNoNumbers[:-1]
    # Let's find what object is it
    for line in open(path.join("references", "list_of_objects.dat")):
        objName = line.strip()
        if fNameNoFilt.startswith(objName):
            addString = fNameNoFilt[len(objName):]
            return objName, filtName, addString
    # No sutable object found
    return None, None, None


class MainApplication(Tk.Frame):
    def __init__(self, *args, **kwargs):
        self.currentObject = None
        self.currentFilter = None
        self.rawImages = []
        self.magData = []
        # cache flats
        self.flats = {}
        for filt in 'bvriXY':
            pathToFlat = path.join("tmp2", "flats", "ff%s.fts" % (filt))  # FIXME
            fData = pyfits.open(pathToFlat)[0].data
            self.flats[filt.lower()] = fData / np.mean(fData)
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
        self.root.after(1000, self.cycle)
        

    def check_out_files(self):
        allFiles = glob.glob(path.join(self.dirName, "*"))
        lightFiles = []
        # get rid of dark files and subdirectories
        for f in allFiles:
            fName = path.basename(f)
            if (not "dark" in fName) and (not "bias" in fName) and (path.isfile(f)):
                lightFiles.append(f)
        newestFile = max(lightFiles, key=path.getctime)
        fNameWithoutPath = path.basename(newestFile)
        # Let's find what object is it
        self.objName, self.filtName, self.addString = parse_object_file_name(fNameWithoutPath)
        if (not self.filtName is None) and (self.filtName.lower() in ("x", "y")):
            self.polarMode = True
        else:
            self.polarMode = False


    def run_computation(self):
        """ This is the main function. It is been
        called as soon as the object is selected."""
        if self.filtName != self.currentFilter:
            # new filter being processed
            self.rawImages = []
            self.darkCleanImages = []
            self.magData = []
            self.currentFilter = self.filtName
        if self.objName != self.currentObject:
            # new object is being processed
            self.rawImages = []
            self.magData = []
            self.darkCleanImages = []
            self.currentObject = self.objName
            self.ref = Reference(self.objName)
            # Clear working directory
            for f in glob.glob(path.join("workDir", "*")):
                if (not "dark" in path.basename(f)) and (not "bias" in path.basename(f)):
                    os.remove(f)

        # Create the list of images to be processed
        newRawImages = []
        for img in sorted(glob.glob(path.join(self.dirName, "%s%s%s*" % (self.objName, self.addString, self.filtName)))):
            if not img in self.rawImages:
                newRawImages.append(img)
        if len(newRawImages) == 0:
            print "no new images"
            return
            
        self.rawImages.extend(newRawImages)

        # Subtract bias and dark files
        self.darkCleanImages.extend(fix_for_bias_dark_and_flat(self.dirName, newRawImages, self.flats[self.filtName.lower()]))
        self.rightPanel.update_message("Bias and dark", "Ok")

        # Coadd images
        if not self.polarMode:
            numOfCoaddedImages = coadd_images(self.darkCleanImages, None)
        else:
            numOfCoaddedImages = coadd_images(self.darkCleanImages, self.filtName)
        self.rightPanel.update_message("Images summed", "%i" % numOfCoaddedImages)

        # Subtract background
        clean_background()
        backCleanFile = path.join("workDir", "back_clean.fits")
        self.imagPanel.plot_fits_file(backCleanFile)

        # Create a catalogue of objects on summed images
        catName = path.join("workDir", "field.cat")
        if self.polarMode:
            # in polar mode we have two catalogues: polar and filtred one.
            catNamePolar = path.join("workDir", "field_polar.cat")
            call_SE(backCleanFile, catNamePolar)
            filterPolCat(catNamePolar, catName, self.filtName)
            self.seCatPolar = SExCatalogue(catNamePolar)
            self.seCat = SExCatalogue(catName)
            # Match objects from reference image on the observed one
            returnCode = self.ref.match_objects(backCleanFile, self.seCatPolar, self.filtName)
        else:
            call_SE(backCleanFile, catName)
            self.seCat = SExCatalogue(catName)
            # Match objects from reference image on the observed one
            returnCode = self.ref.match_objects(backCleanFile, self.seCat)

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
            objSn, objMag, stSn = get_photometry(self.seCat, self.ref, self.filtName)
            self.rightPanel.show_photometry_data(objSn, objMag, stSn)
            self.magData.append((numOfCoaddedImages, objMag))
            # self.graphPanel.plot_magnitude(self.magData)

        elif self.polarMode:
            objSN, objPairSN, stSnList = get_photometry_polar_mode(self.seCatPolar, self.ref)
            self.rightPanel.show_photometry_data_polar_mode(objSN, objPairSN, stSnList)
            # self.rightPanel.remove_objects_from_plot(upDateFigure=True)

        return


[os.remove(f) for f in glob.glob("workDir/*.*")]
MainApplication()




#fix_for_bias_and_dark("./tmpdata/", glob.glob("./tmpdata/s50716*r0*"))


#coadd_images(glob.glob("./tmpdata/s50716v0*"))

# #clean_background()

# ref = Reference("s50716")
# tr = find_shift(ref.fitsFile, "./workDir/summed.fits")
# print tr.v
# print find_objects(ref, tr)

