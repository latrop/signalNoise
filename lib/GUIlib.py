#! /usr/bin/env python

import os
from math import cos, sin, radians
from os import path
import glob
from collections import OrderedDict
import tkinter as Tk
from tkinter import font as tkFont
from tkinter import filedialog as tkFileDialog
from tkinter import messagebox as tkMessageBox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pylab
try:
    import pyfits
except ImportError:
    from astropy.io import fits as pyfits
import numpy as np

from reduction import parse_object_file_name


class MenuBar(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.menubar = Tk.Menu(window.root)
        self.fileMenu = Tk.Menu(self.menubar, tearoff=0)
        self.fileMenu.add_command(label="Select folder", command=self.select_folder)
        self.fileMenu.add_command(label="Rename files", command=self.rename_files)
        self.fileMenu.add_command(label="Reset", command=self.window.reset_new_object)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)

        # background type configuration
        self.backMenu = Tk.Menu(self.menubar, tearoff=0)
        self.backTypeVar = Tk.StringVar()
        self.backTypeVar.set("GLOBAL")
        self.backMenu.add_radiobutton(label="GLOBAL", variable=self.backTypeVar)
        self.backMenu.add_radiobutton(label="LOCAL", variable=self.backTypeVar)
        self.menubar.add_cascade(label="Background", menu=self.backMenu)

        # show menu
        self.showMenu = Tk.Menu(self.menubar, tearoff=0)
        self.showMenu.add_checkbutton(label="Hot pixels", onvalue=True,
                                      offvalue=False, variable=self.window.showHotPixels)
        self.menubar.add_cascade(label="Show", menu=self.showMenu)

        # Alarm menu button
        self.menubar.add_command(label="Alarm", command=self.set_alarm)

        # Polar check menu button
        self.menubar.add_command(label="PolarCheck", command=self.polar_check)

        # Log menu button
        self.logMenu = Tk.Menu(self.menubar, tearoff=0)
        self.logMenu.add_command(label="Show log", command=self.show_log)
        self.logMenu.add_command(label="Select object", command=self.select_object)
        self.menubar.add_cascade(label="History", menu=self.logMenu)
        # Quit
        self.menubar.add_command(label="Quit", command=self.window.on_closing)
        self.window.root.config(menu=self.menubar)

    def select_folder(self):
        # Find newest directory in telescopedata as an initial dir
        if os.name == "nt":
            pathToData = path.join("C:\\", "TelescopeData", "*")
            listOfDirs = filter(path.isdir, glob.glob(pathToData))
            initialdir = max(listOfDirs, key=path.getctime)
        else:
            initialdir = None
        dirPath = tkFileDialog.askdirectory(parent=self.window.root,
                                            title="Open data folder",
                                            initialdir=initialdir)
        dirPath = path.normpath(dirPath)
        self.window.setup(dirPath)

    def set_alarm(self):
        AlarmPopup(self.window)

    def polar_check(self):
        if self.window.filtName.lower() in "xy":
            tkMessageBox.showwarning("PolarCheck",
                                     "Works only in photometry mode.")
            return
        PolarChecker(self.window)

    def rename_files(self):
        RenameFilesPopup(self.window)

    def show_log(self):
        LogWindow(self.window)

    def select_object(self):
        SelectObjectWindow(self.window)


class ImagPanel(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.mainGraph = pylab.Figure(figsize=(6, 4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.mainGraph, master=self.window.root)
        self.fig = self.mainGraph.add_subplot(111)
        self.fig.axes.set_xticks([])
        self.fig.axes.set_yticks([])
        self.fitsPlotInstance = None
        self.objPlotInstance = None
        self.standartsPlotIntance = []
        self.objPairPlotInstance = None
        self.standartsPairPlotInstance = []
        self.hotPixelsPlotInstance = []

        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=0, row=0)

    def remove_objects_from_plot(self, upDateFigure=False):
        if self.objPlotInstance:
            self.objPlotInstance.remove()
            self.objPlotInstance = None
        while self.standartsPlotIntance:
            self.standartsPlotIntance.pop().remove()
        if self.objPairPlotInstance:
            self.objPairPlotInstance.remove()
            self.objPairPlotInstance = None
        while self.standartsPairPlotInstance:
            self.standartsPairPlotInstance.pop().remove()
        while self.hotPixelsPlotInstance:
            self.hotPixelsPlotInstance.pop().remove()
        if upDateFigure:
            self.canvas.show()

    def plot_fits_file(self, fitsName):
        if self.fitsPlotInstance:
            # remove previous plot if exists
            self.fitsPlotInstance.remove()
            self.fitsPlotInstance = None
        hdu = pyfits.open(fitsName)
        data = hdu[0].data.copy()
        hdu.close()
        ySize, xSize = data.shape
        meanValue = np.mean(data)
        stdValue = np.std(data)
        maxValue = meanValue + 2*stdValue
        self.fitsPlotInstance = self.fig.imshow(data, interpolation='gaussian', cmap='gray',
                                                vmin=meanValue, vmax=maxValue)
        self.fig.axis([0, xSize, ySize, 0])
        self.canvas.show()

    def plot_objects(self, reference, polarMode=None, hotPixels=[]):
        self.remove_objects_from_plot()
        if reference.objSEParams is not None:
            xCoord = reference.objSEParams["X_IMAGE"] - 1
            yCoord = reference.objSEParams["Y_IMAGE"] - 1
            self.objPlotInstance = self.fig.plot([xCoord], [yCoord], marker="o", markerfacecolor="none",
                                                 markersize=15, markeredgewidth=2, markeredgecolor="r")[0]
        else:
            self.objPlotInstance = self.fig.plot([reference.xObjObs-1], [reference.yObjObs-1], marker="o",
                                                 markerfacecolor="none", markersize=15, markeredgewidth=2,
                                                 markeredgecolor="0.75")[0]

        for st in reference.standartsObs:
            if st['seParams'] is None:
                stx = st['xCen']-1
                sty = st['yCen']-1
                markColor = "0.75"
            else:
                stx = st['seParams']["X_IMAGE"]-1
                sty = st['seParams']["Y_IMAGE"]-1
                markColor = "g"
            self.standartsPlotIntance.append(self.fig.plot([stx], [sty], marker="o", markerfacecolor="none",
                                                           linestyle="", markersize=15, markeredgewidth=2,
                                                           markeredgecolor=markColor)[0])

        if polarMode:
            if reference.objPairSEParams is not None:
                self.objPairPlotInstance = self.fig.plot([reference.objPairSEParams["X_IMAGE"]-1],
                                                         [reference.objPairSEParams["Y_IMAGE"]-1],
                                                         marker="o", markerfacecolor="none", markersize=15,
                                                         markeredgewidth=2, markeredgecolor="r")[0]
            else:
                self.objPairPlotInstance = self.fig.plot([reference.xObjPairObs-1],
                                                         [reference.yObjPairObs-1],
                                                         marker="o", markerfacecolor="none", markersize=15,
                                                         markeredgewidth=2, markeredgecolor="0.75")[0]

            for st in reference.standartPairsObs:
                if st['seParams'] is None:
                    stx = st['xCen']-1
                    sty = st['yCen']-1
                    markColor = "0.75"
                else:
                    stx = st['seParams']["X_IMAGE"]-1
                    sty = st['seParams']["Y_IMAGE"]-1
                    markColor = "g"
                self.standartsPairPlotInstance.append(self.fig.plot([stx], [sty], marker="o", markerfacecolor="none",
                                                                    linestyle="", markersize=15, markeredgewidth=2,
                                                                    markeredgecolor="g")[0])

        # plot hot pixels
        if self.window.showHotPixels.get():
            self.hotPixelsPlotInstance.append(self.fig.plot(hotPixels[1], hotPixels[0], marker="x", markersize=5,
                                                            markeredgecolor="b", markerfacecolor="b",
                                                            markeredgewidth=1, linestyle="")[0])
        self.canvas.show()


class RightPanel(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.panel = Tk.Frame(self.window.root)
        self.panel.grid(column=1, row=0)
        self.objectInfoLabelValue = Tk.StringVar()
        Tk.Label(self.panel, textvariable=self.objectInfoLabelValue).grid(column=0, row=0, sticky=Tk.W)
        self.objectInfoLabelValue.set("Object is not selected.")
        self.messagesString = Tk.StringVar()
        self.messages = {}
        Tk.Label(self.panel, textvariable=self.messagesString).grid(column=0, row=2, sticky=Tk.W)
        self.photometryString = Tk.StringVar()
        Tk.Label(self.panel, textvariable=self.photometryString).grid(column=0, row=3, sticky=Tk.W)
        self.badImagesString = Tk.StringVar()
        Tk.Label(self.panel, textvariable=self.badImagesString, fg="red").grid(column=0, row=4, sticky=Tk.W)

    def update_object_info(self, objName, filtName, addString):
        s = "Object: %s" % (objName)
        if addString:
            s += "(%s)" % (addString)
        s += ", filter %s" % (filtName)
        self.objectInfoLabelValue.set(s)

    def update_message(self, key, msg):
        if (key in self.messages) and (self.messages[key] == msg):
            # nothing new
            return False
        self.messages[key] = msg
        self.messagesString.set("\n"+"\n".join(["%s: %s" % (key, self.messages[key])
                                                for key in self.messages if self.messages[key]]))
        if "Error" in self.messages:
            self.photometryString.set("")
        return True

    def show_photometry_data(self, objSn, objMag, objMagSigma, stSn):
        string = ""
        if objMag is None:
            string += "Mag: undef\n"
        else:
            string += u"Mag: %1.2f \u00B1 %1.2f\n" % (objMag, objMagSigma)
        if objSn is None:
            string += "Obj S/N: undef\n"
        else:
            string += "Obj S/N: %1.0f\n" % (objSn)
        for key in sorted(stSn):
            value = stSn[key]
            if value is None:
                string += "%s S/N: undef\n" % (key)
            else:
                string += "%s S/N: %1.0f\n" % (key, value)
        self.photometryString.set(string)

    def show_photometry_data_polar_mode(self, objSN, objPairSN, stSnList, fluxRatios):
        string = "Obj S/N: "
        if objSN is None:
            string += "undef/"
        else:
            string += "%1.0f/" % (objSN)
        if objPairSN is None:
            string += "undef\n"
        else:
            string += "%1.0f\n" % (objPairSN)

        for key in sorted(stSnList):
            values = stSnList[key]
            string += "%s S/N: " % (key)
            if values[0] is None:
                string += "undef/"
            else:
                string += "%1.0f/" % (values[0])
            if values[1] is None:
                string += "undef\n"
            else:
                string += "%1.0f\n" % (values[1])

        string += "\n"
        if fluxRatios["obj"] is not None:
            string += "Obj FR: %1.3f\n" % fluxRatios["obj"]
        else:
            string += "Obj FR: Undef\n"
        if fluxRatios["st"] is not None:
            string += "St FR: %1.3f\n" % fluxRatios["st"]
        else:
            string += "St FR: undef\n"

        self.photometryString.set(string)

    def update_bad_images_info(self, badImagesList):
        if badImagesList:
            string = "\nBad images:\n"
            string += "\n".join(badImagesList)
        else:
            string = ""
        self.badImagesString.set(string)


class AlarmPopup(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.top = Tk.Toplevel(window.root)
        self.top.geometry('+%i+%i' % (window.root.winfo_x()+30, window.root.winfo_y()+30))
        Tk.Label(self.top, text="Exposures").grid(column=0, row=0, padx=5, pady=5)
        self.entry = Tk.Entry(self.top, width=5)
        self.entry.insert(0, str(abs(window.desiredExposures)))
        self.entry.grid(column=1, row=0, padx=5, pady=5)
        self.button = Tk.Button(self.top, text="OK", command=self.ok)
        self.button.grid(column=0, row=1, padx=5, pady=5)
        self.cancelButton = Tk.Button(self.top, text="Cancel", command=self.top.destroy)
        self.cancelButton.grid(column=1, row=1, padx=5, pady=5)

    def ok(self):
        self.window.desiredExposures = int(self.entry.get())
        self.top.destroy()


class RenameFilesPopup(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.top = Tk.Toplevel(window.root)
        self.top.geometry('+%i+%i' % (window.root.winfo_x()+30, window.root.winfo_y()+30))
        Tk.Label(self.top, text="Input desired number of exposures").grid(column=0, row=0, padx=5, pady=5)
        self.entry = Tk.Entry(self.top, width=5)
        self.entry.insert(0, str(abs(window.numOfCoaddedImages)))
        self.entry.grid(column=1, row=0, padx=5, pady=5)
        self.button = Tk.Button(self.top, text="OK", command=self.ok)
        self.button.grid(column=0, row=1, padx=5, pady=5)
        self.cancelButton = Tk.Button(self.top, text="Cancel", command=self.top.destroy)
        self.cancelButton.grid(column=1, row=1, padx=5, pady=5)
        self.top.wm_attributes("-topmost", 1)

    def ok(self):
        desiredExposures = int(self.entry.get())
        if desiredExposures >= len(self.window.rawImages):
            minorExp = self.window.rename_files(desiredExposures)
            self.window.reset_new_object()
            tkMessageBox.showinfo("Rename files",
                                  "All done, minor exposure is %i" % (minorExp+1))
            self.top.destroy()
        else:
            errorString = "The number of desired exposures should be bigger than the number of already taken ones."
            tkMessageBox.showwarning("Rename files", errorString)


class PolarChecker(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.top = Tk.Toplevel(window.root)
        # self.top.geometry("500x1000")
        self.yFitsPlotInstance = None
        self.xFitsPlotInstance = None
        self.yObjPlotInstance = None
        self.xObjPlotInstance = None
        self.yStdPlotInstance = []
        self.xStdPlotInstance = []
        self.rotation = 0.0
        self.cosa = 1
        self.sina = 0
        # Initialisation of graph for y-mode
        self.yGraph = pylab.Figure(figsize=(6, 4), dpi=80)
        self.yCanvas = FigureCanvasTkAgg(self.yGraph, master=self.top)
        self.yFig = self.yGraph.add_subplot(111)
        self.yFig.axes.set_xticks([])
        self.yFig.axes.set_yticks([])
        self.yFig.axis([0, 382, 255, 0])
        self.yFig.axes.set_title("Y-mode")
        self.yCanvas.show()
        self.yCanvas.get_tk_widget().grid(column=0, row=0, columnspan=4)

        # Initialisation of graph for x-mode
        self.xGraph = pylab.Figure(figsize=(6, 4), dpi=80)
        self.xCanvas = FigureCanvasTkAgg(self.xGraph, master=self.top)
        self.xFig = self.xGraph.add_subplot(111)
        self.xFig.axes.set_xticks([])
        self.xFig.axes.set_yticks([])
        self.xFig.axis([0, 382, 255, 0])
        self.xFig.axes.set_title("X-mode")
        self.xCanvas.show()
        self.xCanvas.get_tk_widget().grid(column=0, row=1, columnspan=4)

        # Add some bottons
        self.rotIncButton = Tk.Button(self.top, text="+5 deg", command=lambda: self.rotate(5.0))
        self.rotIncButton.grid(column=0, row=2, pady=10)
        self.rotDecButton = Tk.Button(self.top, text="-5 deg", command=lambda: self.rotate(-5.0))
        self.rotDecButton.grid(column=2, row=2, pady=10)
        self.okButton = Tk.Button(self.top, text="Ok", command=self.top.destroy)
        self.okButton.grid(column=3, row=2, pady=10)

        # Text with current rotation
        self.angleValue = Tk.StringVar()
        self.angleValue.set("0")
        Tk.Label(self.top, textvariable=self.angleValue).grid(column=1, row=2)

        # show catalogue
        self.show_cat()

    def rotate(self, angle):
        self.rotation += angle
        self.angleValue.set("%i" % (int(self.rotation)))
        self.cosa = cos(radians(self.rotation))
        self.sina = sin(radians(self.rotation))
        # Show catalogues
        self.show_cat()

    def clear_fig(self):
        # remove previous plots if exist
        if self.yFitsPlotInstance is not None:
            self.yFitsPlotInstance.remove()
            self.yFitsPlotInstance = None

        if self.xFitsPlotInstance is not None:
            self.xFitsPlotInstance.remove()
            self.xFitsPlotInstance = None

        if self.yObjPlotInstance is not None:
            self.yObjPlotInstance.remove()
            self.yObjPlotInstance = None

        if self.xObjPlotInstance is not None:
            self.xObjPlotInstance.remove()
            self.xObjPlotInstance = None

        while self.yStdPlotInstance:
            self.yStdPlotInstance.pop().remove()
        while self.xStdPlotInstance:
            self.xStdPlotInstance.pop().remove()

        # self.yCanvas.show()
        # self.xCanvas.show()

    def show_cat(self):
        """ Show catalogue as an image to check if some polar pairs are too close"""
        self.clear_fig()
        xSize = 382
        ySize = 255
        xCen = 191
        yCen = 127.5

        gridX, gridY = np.meshgrid(range(xSize), range(ySize))

        # Create images for x and y mode
        dataY = np.zeros((ySize, xSize))
        dataX = np.zeros((ySize, xSize))
        for obj in self.window.seCat:
            if (obj["FLUX_AUTO"] <= 0) or (obj["FWHM_IMAGE"] <= 0):
                continue
            xObjOrig = obj["X_IMAGE"]
            yObjOrig = obj["Y_IMAGE"]
            # Lets rotate points by angle self.rotation around the centre
            xObj = self.cosa * (xObjOrig-xCen) - self.sina * (yObjOrig-yCen) + xCen
            yObj = self.sina * (xObjOrig-xCen) + self.cosa * (yObjOrig-yCen) + yCen

            xPairY = xObj + 17.7
            yPairY = yObj - 0.7
            xPairX = xObj + 12.5
            yPairX = yObj - 12.5

            r = 3*int(obj["FWHM_IMAGE"])
            # set indexes to work only for part of the image
            idxObj = np.s_[int(yObj)-r:int(yObj)+r+1,
                           int(xObj)-r:int(xObj)+r+1]
            idxPairY = np.s_[int(yPairY)-r:int(yPairY)+r+1,
                             int(xPairY)-r:int(xPairY)+r+1]
            idxPairX = np.s_[int(yPairX)-r:int(yPairX)+r+1,
                             int(xPairX)-r:int(xPairX)+r+1]

            sqDistsObj = np.zeros((ySize, xSize))
            sqDistsPairY = np.zeros((ySize, xSize))
            sqDistsPairX = np.zeros((ySize, xSize))

            sqDistsObj[idxObj] = (gridX[idxObj]-xObj)**2.0 + (gridY[idxObj]-yObj)**2.0
            sqDistsPairY[idxPairY] = (gridX[idxPairY]-xPairY)**2.0 + (gridY[idxPairY]-yPairY)**2.0
            sqDistsPairX[idxPairX] = (gridX[idxPairX]-xPairX)**2.0 + (gridY[idxPairX]-yPairX)**2.0

            objImag = np.zeros((ySize, xSize))
            objImag[idxObj] = obj["FLUX_AUTO"] * np.exp(-sqDistsObj[idxObj]/(2*obj["FWHM_IMAGE"]))
            dataY[idxObj] += objImag[idxObj]
            dataY[idxPairY] += obj["FLUX_AUTO"] * np.exp(-sqDistsPairY[idxPairY]/(2*obj["FWHM_IMAGE"]))
            dataX[idxObj] += objImag[idxObj]
            dataX[idxPairX] += obj["FLUX_AUTO"] * np.exp(-sqDistsPairX[idxPairX]/(2*obj["FWHM_IMAGE"]))

        yMean = np.mean(dataY)
        yStd = np.std(dataY)
        self.yFitsPlotInstance = self.yFig.imshow(dataY, interpolation='gaussian',
                                                  cmap='gray', vmin=yMean, vmax=yMean+2*yStd)

        xMean = np.mean(dataX)
        xStd = np.std(dataX)
        self.xFitsPlotInstance = self.xFig.imshow(dataX, interpolation='gaussian',
                                                  cmap='gray', vmin=xMean, vmax=xMean+2*xStd)

        # Overplot location of the object and reference stars
        if self.window.ref.objSEParams is not None:
            objX0 = self.window.ref.objSEParams["X_IMAGE"]
            objY0 = self.window.ref.objSEParams["Y_IMAGE"]
        else:
            objX0 = self.window.ref.xObjObs
            objY0 = self.window.ref.yObjObs
        objX = self.cosa * (objX0-xCen) - self.sina * (objY0-yCen) + xCen
        objY = self.sina * (objX0-xCen) + self.cosa * (objY0-yCen) + yCen
        self.yObjPlotInstance = self.yFig.plot([objX, objX+17.7], [objY, objY-0.7], marker="o", markerfacecolor="none",
                                               markersize=15, markeredgewidth=2, markeredgecolor="r", linestyle="")[0]
        self.xObjPlotInstance = self.xFig.plot([objX, objX+12.5], [objY, objY-12.5], marker="o", markerfacecolor="none",
                                               markersize=15, markeredgewidth=2, markeredgecolor="r", linestyle="")[0]

        for st in self.window.ref.standartsObs:
            if st['seParams'] is not None:
                stx0 = st['seParams']["X_IMAGE"]-1
                sty0 = st['seParams']["Y_IMAGE"]-1
                markColor = "g"
            else:
                stx0 = st['xCen'] - 1
                sty0 = st['yCen'] - 1
                markColor = "0.75"
            stx = self.cosa * (stx0-xCen) - self.sina * (sty0-yCen) + xCen
            sty = self.sina * (stx0-xCen) + self.cosa * (sty0-yCen) + yCen
            self.yStdPlotInstance.append(self.yFig.plot([stx, stx+17.7], [sty, sty-0.7], marker="o",
                                                        markerfacecolor="none", linestyle="", markersize=15,
                                                        markeredgewidth=2, markeredgecolor=markColor)[0])
            self.xStdPlotInstance.append(self.xFig.plot([stx, stx+12.5], [sty, sty-12.5], marker="o",
                                                        markerfacecolor="none", linestyle="", markersize=15,
                                                        markeredgewidth=2, markeredgecolor=markColor)[0])

        self.yCanvas.show()
        self.xCanvas.show()


class LogWindow(Tk.Frame):
    """
    Window to show photometry log data
    """
    def __init__(self, window):
        self.top = Tk.Toplevel(window.root)
        self.textFrame = Tk.Text(self.top)
        self.scroll = Tk.Scrollbar(self.top)
        self.scroll.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.textFrame.pack(side=Tk.LEFT, fill=Tk.Y)
        self.scroll.config(command=self.textFrame.yview)
        self.textFrame.config(yscrollcommand=self.scroll.set,
                              font=tkFont.Font(family="Times", size=12))
        # Put text into frame
        for key in window.photoLog:
            objName = key.split(":")[0]
            addString = key.split(":")[1]
            if addString:
                logString = "%s(%s): " % (objName, addString)
            else:
                logString = "%s:  " % (objName)
            for filtName in "bvri":
                if filtName in window.photoLog[key]:
                    magValue = window.photoLog[key][filtName]
                    if magValue is not None:
                        logString += "M_%s=%1.2f  " % (filtName, magValue)
            logString += "\n"
            self.textFrame.insert(Tk.END, logString)


class SelectObjectWindow(Tk.Frame):
    """
    Window to select previously observed objects to be processed
    """
    def __init__(self, window):
        self.window = window
        self.top = Tk.Toplevel(window.root)
        self.top.protocol('WM_DELETE_WINDOW', self.close)

        # Before going into processing let's check the directory real quick to
        # find the list of already observed objects.
        self.list_of_all_observed_objects = OrderedDict()
        allFiles = glob.glob(path.join(self.window.dirName, "*.FIT"))
        allFiles.sort(key=path.getctime)
        for fName in allFiles:
            fNameWithoutPath = path.basename(fName)
            if ("dark" not in fNameWithoutPath) and ("bias" not in fNameWithoutPath) and (path.isfile(fName)):
                objName, filtName, addString, frameNumber = parse_object_file_name(fNameWithoutPath)
                objStr = "%s:%s" % (objName, addString)
                if objStr not in self.list_of_all_observed_objects:
                    # A new object found -> add it to the dictionary
                    self.list_of_all_observed_objects[objStr] = []
                if filtName not in self.list_of_all_observed_objects[objStr]:
                    # Add a new filter for the filterlist of the current object
                    self.list_of_all_observed_objects[objStr].append(filtName)

        # Scrollbar
        self.scrollbar = Tk.Scrollbar(self.top)
        self.scrollbar.grid(column=1, row=0, rowspan=8, sticky=Tk.S+Tk.N)

        # List of all observed objects
        self.objectsListBox = Tk.Listbox(self.top, selectmode=Tk.SINGLE,
                                         height=10, yscrollcommand=self.scrollbar.set)
        self.objectsListBox.grid(column=0, row=0, rowspan=8)
        for objStr in self.list_of_all_observed_objects.keys():
            self.objectsListBox.insert(Tk.END, objStr)
        self.objectsListBox.selection_set(Tk.END)
        self.objectsListBox.bind('<<ListboxSelect>>', self.select_object)
        self.scrollbar.config(command=self.objectsListBox.yview)
        self.objectsListBox.yview_moveto(1)

        # Radoibuttons with filtres
        self.selectedFilter = Tk.StringVar()
        self.selectedFilter.set("r")
        self.selectedFilter.trace("w", self.compute_previous_object)
        filters = ["b", "v", "r", "i", "y", "x"]
        self.filterRadioButtons = {}
        for i, filt in enumerate(filters):
            self.filterRadioButtons[filt] = Tk.Radiobutton(self.top, text=filt, variable=self.selectedFilter,
                                                           value=filt)
            self.filterRadioButtons[filt].grid(column=2, row=i)

        # Close button
        self.closeButton = Tk.Button(self.top, text="Close", command=self.close)
        self.closeButton.grid(column=2, row=len(filters)+1)
        self.objectsListBox.event_generate("<<ListboxSelect>>")

    def select_object(self, event):
        self.objStr = self.objectsListBox.get(self.objectsListBox.curselection())
        observed_filters = self.list_of_all_observed_objects[self.objStr]
        # Make visible only radiobuttons that correspond to observed filters
        filters = ["b", "v", "r", "i", "y", "x"]
        for filt in filters:
            if filt in observed_filters:
                self.filterRadioButtons[filt].config(state="normal")
            else:
                self.filterRadioButtons[filt].config(state="disabled")
        # Select r filter if it was observed, otherwise select just first observed filter
        if 'r' in observed_filters:
            self.selectedFilter.set('r')
        else:
            self.selectedFilter.set(observed_filters[0])

    def compute_previous_object(self, *args):
        objName, addString = self.objStr.split(":")
        self.window.object_selected_manually = True
        self.window.objName = objName
        self.window.filtName = self.selectedFilter.get()
        self.window.addString = addString

    def close(self):
        self.window.object_selected_manually = False
        self.top.destroy()
