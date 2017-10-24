#! /usr/bin/env python

import glob
from os import path
from os import remove

import numpy
try:
    import pyfits
except ImportError:
    from astropy.io import fits as pyfits

import alipylocal as alipy


class Reference(object):
    """ Class stores reference data (image and
    coordinates of object and standart stars
    within this image"""
    def __init__(self, objName):
        """
        There can be more than one reference file for a given object.
        So the strategy is to try different references untill good match
        obtained and then use correspondent coords.dat file.
        """
        self.refImages = sorted(glob.glob(path.join("references", objName, "ref*.cat")))
        self.coordMagDataFiles = sorted(glob.glob(path.join("references", objName, "coords*.dat")))
        self.transform = None
        self.standarts = []

        # Check if there are aperture settings
        aperPath = path.join("references", objName, "aperture.dat")
        if path.exists(aperPath):
            aperFile = open(aperPath)
            aperData = aperFile.readlines()[0]
            self.apertureSize = float(aperData.split()[1])
            aperFile.close()
            print "+++++++ aperture = %1.1f" % self.apertureSize
        else:
            self.apertureSize = None

        
    def load_coord_mags(self, nRef):
        """
        Function loads reference coordinates and magnitudes from
        coords.dat file that corresponds to a successful reference.
        nRef is the ordinal number of this successful reference.
        """
        for line in open(self.coordMagDataFiles[nRef]):
            if line.strip().startswith("#"):
                continue
            name = line.split()[0]
            xCen = float(line.split()[1])
            yCen = float(line.split()[2])
            if name == "obj":  # coordinates of object
                self.xObj = xCen
                self.yObj = yCen
                self.objSEParams = None
            elif name.startswith("st"):
                self.standarts.append({"name": name, "xCen": xCen, "yCen": yCen})
                for i, filt in enumerate("bvri"):
                    # load magnitudes of standarts
                    magStr = line.split()[i+3]
                    magVal = float(magStr) if "-" not in magStr else None
                    self.standarts[-1]["mag%s" % filt] = magVal

    def find_shift(self, observedImage, polarMode):
        if not path.exists(observedImage):
            print "%s not found!" % (observedImage)
            return
        for nRef, refImage in enumerate(self.refImages):
            # Here we try different references to find the
            # one that matches
            ident = alipy.ident.run(refImage, [observedImage], visu=False, verbose=True, polarMode=polarMode)
            if (ident is not None) and (ident[0].ok is True):
                # Good reference is found
                self.transform = ident[0].trans.inverse()
                if not self.standarts:
                    self.load_coord_mags(nRef)
                return
        # If we run out of reference files, but proper match was not
        # found, then we can use the last good transformation (obtained
        # from prevous file set of in different filter). If there was
        # not previous good transformation, then the algorithm
        # will use default (None) transformation and do nothing
        print "not ok"
            
    def apply_transform(self):
        """ Function finds coordinates of object and standarts
        on the observed image"""
        self.xObjObs, self.yObjObs = self.transform.apply((self.xObj, self.yObj))

        # find coordinates of standarts
        self.standartsObs = []
        for st in self.standarts:
            x, y = self.transform.apply((st["xCen"], st["yCen"]))
            self.standartsObs.append({"name": st['name'], "xCen": x, "yCen": y})            
            
    def match_objects(self, observedImage, catalogue, polarMode=None, matchOnly=False):
        # Let's find where on the observed image objects should
        # be according to reference image using affine transform
        if not matchOnly:
            self.find_shift(observedImage, polarMode=polarMode)
            if self.transform is None:
                # No transform found
                return 1
            self.apply_transform()

        # Now we need to find for objects their sextractor-generated parameters
        self.objSEParams = catalogue.find_nearest(self.xObjObs, self.yObjObs)
        for st in self.standartsObs:
            st["seParams"] = catalogue.find_nearest(st['xCen'], st['yCen'])

        # If polar mode, we have to find polar counterparts for every object
        if polarMode is None:
            self.xObjPairObs = None
            self.yObjPairObs = None
            self.standartPairsObs = None
        else:
            self.find_polar_pairs(polarMode)
            self.objPairSEParams = catalogue.find_nearest(self.xObjPairObs, self.yObjPairObs)
            for st in self.standartPairsObs:
                st["seParams"] = catalogue.find_nearest(st['xCen'], st['yCen'])
        return 0


    def find_polar_pairs(self, polarMode):
        self.standartPairsObs = []
        if polarMode.lower() == 'x':
            self.xObjPairObs = self.xObjObs + 12.5
            self.yObjPairObs = self.yObjObs - 12.5
            for st in self.standartsObs:
                x = st["xCen"] + 12.5
                y = st["yCen"] - 12.5
                self.standartPairsObs.append({"name": "%s_b" % (st['name']),
                                              "xCen": x, "yCen": y})
        elif polarMode.lower() == 'y':
            self.xObjPairObs = self.xObjObs + 17.7
            self.yObjPairObs = self.yObjObs - 0.7
            for st in self.standartsObs:
                x = st["xCen"] + 17.7
                y = st["yCen"] - 0.7
                self.standartPairsObs.append({"name": "%s_b" % (st['name']),
                                              "xCen": x, "yCen": y})

    def get_standatds_fwhm(self):
        fwhms = [st["seParams"]["FWHM_IMAGE"] for st in self.standartsObs if st["seParams"] is not None]
        return numpy.mean(fwhms)


def coadd_images(imageList, polarMode):
    """ Function shifts images from imageList, so they
    all match, and then coadds them to one """
    # First step: match images
    refImage = imageList[0]
    outputshape = alipy.align.shape(refImage)
    # We don't want to process images that was already shifted, so we
    # should check for such images in our list. This can be done by
    # searching for "*affineremap" images in workDir directory
    imagesToAlign = []
    imagesToCoadd = []
    for image in imageList[1:]:
        imgName = path.splitext(path.basename(image))[0]
        pathToFile = path.join("workDir", "%s_affineremap.fits" % (imgName))
        if not path.exists(pathToFile):
            imagesToAlign.append(image)
        else:
            # if there is already an aligned image, then we can just
            # add it to imagesToCoadd list
            imagesToCoadd.append(pathToFile)
    # Shifting should be done only if there is new images:
    if imagesToAlign:
        print ("%s" % str(polarMode)) * 10
        identifications = alipy.ident.run(refImage, imagesToAlign, visu=False, verbose=True, r=10.0,
                                          polarMode=polarMode, refpolar=True)
        for ident in identifications:
            if ident.ok is True:
                # if alignment is ok, then transform image
                alipy.align.affineremap(ident.ukn.filepath, ident.trans, shape=outputshape,
                                        makepng=False, outdir='workDir', verbose=False)
                # and add it to coadd list
                imgName = path.splitext(path.basename(ident.ukn.filepath))[0]
                pathToFile = path.join("workDir", "%s_affineremap.fits" % (imgName))
                imagesToCoadd.append(pathToFile)
            else:
                imagesToCoadd.append(ident.ukn.filepath)
                print "not ok: %s" % (path.basename(ident.ukn.filepath))

    # refImage image was not remapped, since all other images was remapped
    # to match it, but we want to coadd it as well:
    imagesToCoadd.append(refImage)


    # Second step: coadd images    
    data = numpy.zeros((outputshape[1], outputshape[0]))
    maskData = numpy.ones((outputshape[1], outputshape[0]), dtype=float)
    for img in imagesToCoadd:
        hdu = pyfits.open(img)
        data += hdu[0].data
        inds = numpy.where(hdu[0].data == 0.0)
        maskData[inds] = 0.0
        hdu.close()
    # fix unoverlapped regions by setting median value there
    data[numpy.where(maskData==0.0)] = numpy.median(data)
    # save summed data to summed.fits file
    outHdu = pyfits.PrimaryHDU(data=data)
    pathToFile = path.join("workDir", "summed.fits")
    outHdu.writeto(pathToFile, clobber=True)
    return len(imagesToCoadd)

