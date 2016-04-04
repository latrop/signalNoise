#! /usr/bin/env python

import glob
from os import path
from os import remove

import numpy
import pyfits
import alipylocal as alipy


class Reference(object):
    """ Class stores reference data (image and
    coordinates of object and standart stars
    within this image"""
    def __init__(self, objName):
        self.refImage = path.join("references", objName, "ref.cat")
        dataFile = path.join("references", objName, "coords.dat")
        self.standarts = []
        self.transform = None
        if (path.exists(self.refImage)) and (path.exists(dataFile)):
            for line in open(dataFile):
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
                    magb = line.split()[3]
                    if "-" not in magb:
                        magb = float(magb)
                    else:
                        magb = None
                    magv = line.split()[4]
                    if "-" not in magv:
                        magv = float(magv)
                    else:
                        magv = None
                    magr = line.split()[5]
                    if "-" not in magr:
                        magr = float(magr)
                    else:
                        magr = None
                    magi = line.split()[6]
                    if "-" not in magi:
                        magi = float(magi)
                    else:
                        magi = None

                    self.standarts.append({"name": name, "xCen": xCen, "yCen": yCen, "magb": magb,
                                           "magv": magv, "magr": magr, "magi": magi})

    def apply_transform(self):
        """ Function finds coordinates of object and standarts
        on the observed image"""
        matTransform = self.transform.matrixform()
        scaleRotMatrix = matTransform[0]
        translationVector = matTransform[1]
        # x' = x*A + b:
        # find coordinates of the object
        xi = self.xObj-translationVector[0]
        yi = self.yObj-translationVector[1]
        self.xObjObs, self.yObjObs = numpy.dot((xi, yi), scaleRotMatrix)

        # find coordinates of standarts
        self.standartsObs = []
        for st in self.standarts:
            xi = st["xCen"]-translationVector[0]
            yi = st["yCen"]-translationVector[1]
            x, y = numpy.dot((xi, yi), scaleRotMatrix)
            self.standartsObs.append({"name": st['name'], "xCen": x, "yCen": y})


    def find_shift(self, observedImage, polarMode):
        if not path.exists(observedImage):
            print "%s not found!" % (observedImage)
            return
        if not path.exists(self.refImage):
            print "%s not found!" % (self.refImage)
            return
        ident = alipy.ident.run(self.refImage, [observedImage], visu=False, verbose=True, polarMode=polarMode)
        if ident is None:
            print "not ok"
            self.transform = None
        if ident[0].ok is True:
            self.transform = ident[0].trans
        else:
            print "not ok"
            self.transform = None

    def match_objects(self, observedImage, catalogue, polarMode=None):
        # Let's find where on the observed image objects should
        # be according to reference image using affine transform
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
    if len(imagesToAlign) > 0:
        print ("%s" % str(polarMode)) * 10
        identifications = alipy.ident.run(refImage, imagesToAlign, visu=False, verbose=True, r=10.0, polarMode=polarMode, refpolar=True)
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
    # to match it, but we want to coadd in as well:
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
    if path.exists(pathToFile):
        remove(pathToFile)
    outHdu.writeto(pathToFile)
    # maskHDU = pyfits.PrimaryHDU(data=maskData)
    # pathToFile = path.join("workDir", "bad_pixels.fits")
    # if path.exists(pathToFile):
    #     remove(pathToFile)
    # maskHDU.writeto(pathToFile)
    return len(imagesToCoadd)

