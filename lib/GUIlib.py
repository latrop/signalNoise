#! /usr/bin/env python

import os
import itertools
from os import path
import glob
import Tkinter as Tk
import tkFileDialog, tkMessageBox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import pylab
import pyfits
import numpy as np
from scipy.spatial import cKDTree

from lib.alignment import *


class MenuBar(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.menubar = Tk.Menu(window.root)
        self.fileMenu = Tk.Menu(self.menubar, tearoff=0)
        # self.fileMenu.add_command(label="Select object", command=self.select_object)
        self.fileMenu.add_command(label="Select folder", command=self.select_folder)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)

        # background type configuration
        self.backMenu = Tk.Menu(self.menubar, tearoff=0)
        self.backTypeVar = Tk.StringVar()
        self.backTypeVar.set("GLOBAL")
        self.backMenu.add_radiobutton(label="GLOBAL", variable=self.backTypeVar)
        self.backMenu.add_radiobutton(label="LOCAL", variable=self.backTypeVar)
        self.menubar.add_cascade(label="Background", menu=self.backMenu)

        self.menubar.add_command(label="Alarm", command=self.set_alarm)

        self.menubar.add_command(label="PolarCheck", command=self.polar_check)

        self.menubar.add_command(label="Quit", command=self.window.on_closing)
        self.window.root.config(menu=self.menubar)

    # def select_object(self):
    #     fileName = tkFileDialog.askopenfilename(parent=self.window.root,
    #                                             filetypes=[("FITS files", "*.FIT")],
    #                                             title="Open file to load parameters")
    #     # For some reason tkFileDialog returns path in POSIX
    #     # format even in NT system, so we have to fix it
    #     fileName = path.normpath(fileName)
    #     dirPath = path.dirname(fileName)
    #     fNameWithoutPath = path.basename(fileName)
    #     objName, filtName, addString = parse_object_file_name(fNameWithoutPath)
    #     self.window.rightPanel.update_object_info(objName, filtName, addString)
    #     self.window.run_computation(dirPath, objName, filtName, addString)

    def select_folder(self):
        # Find newest directory in telescopedata as an initial dir
        if os.name == "nt":
            pathToData = path.join("C:\\", "TelescopeData", "*")
            listOfDirs = filter(path.isdir, glob.glob(pathToData))
            initialdir  = max(listOfDirs, key=path.getctime)
        else:
            initialdir = None
        dirPath = tkFileDialog.askdirectory(parent=self.window.root,
                                            title="Open data folder",
                                            initialdir=initialdir)
        dirPath = path.normpath(dirPath)
        self.window.setup(dirPath)

    def set_alarm(self):
        popup = AlarmPopup(self.window)

    def polar_check(self):
        PolarChecker(self.window)


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

        # self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.window.root)
        # self.toolbar.update()
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=0, row=0) #pack(side=Tk.TOP)

    def remove_objects_from_plot(self, upDateFigure=False):
        if self.objPlotInstance:
            self.objPlotInstance.remove()
            self.objPlotInstance = None
        if len(self.standartsPlotIntance) > 0:
            for plotInstance in self.standartsPlotIntance:
                plotInstance.remove()
            self.standartsPlotIntance = []
        if self.objPairPlotInstance:
            self.objPairPlotInstance.remove()
            self.objPairPlotInstance = None
        if len(self.standartsPairPlotInstance) > 0:
            for plotInstance in self.standartsPairPlotInstance:
                plotInstance.remove()
            self.standartsPairPlotInstance = []
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
        self.fitsPlotInstance = self.fig.imshow(data, interpolation='nearest', cmap='gray',
                                                vmin=meanValue, vmax=meanValue+2*stdValue)
        self.fig.axis([0, xSize, ySize, 0])
        self.canvas.show()
        #self.canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

    def plot_objects(self, reference, polarMode=None):
        self.remove_objects_from_plot()
        if not reference.objSEParams is None:
            self.objPlotInstance = self.fig.plot([reference.objSEParams["X_IMAGE"]-1], [reference.objSEParams["Y_IMAGE"]-1],
                                                 marker="o", markerfacecolor="none", markersize=15, markeredgewidth=2,
                                                 markeredgecolor="r")[0]
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
            if not reference.objPairSEParams is None:
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
        self.canvas.show()


class GraphPanel(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.mainGraph = pylab.Figure(figsize=(6, 2), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.mainGraph, master=self.window.root)
        self.fig = self.mainGraph.add_subplot(111)
        self.mainGraph.tight_layout()
        self.graphInstance = None
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=0, row=1) #pack(side=Tk.BOTTOM)

    def remove_objects_from_plot(self, upDateFigure=False):
        if self.graphInstance:
            self.graphInstance.remove()
            self.graphInstance = None
            if upDateFigure:
                self.canvas.show()
            
    def plot_magnitude(self, magData):
        self.remove_objects_from_plot()
        x = [m[0] for m in magData]
        y = [m[1] for m in magData]
        mVal = np.mean(y)
        self.fig.axis([x[0]-0.5, x[-1]+0.5, mVal-0.1, mVal+0.1])
        self.graphInstance = self.fig.plot(x, y, "bo", linestyle="-", color="b")[0]
        self.fig.axes.set_xticks([int(i) for i in x])
        self.mainGraph.tight_layout()
        self.canvas.show()


class RightPanel(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.panel = Tk.Frame(self.window.root)
        self.panel.grid(column=1, row=0)# pack(side=Tk.RIGHT, expand=1)
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
        if len(badImagesList) > 0:
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


class PolarChecker(Tk.Frame):
    def __init__(self, window):
        self.rotation = 0.0
        self.window = window
        self.top = Tk.Toplevel(window.root)
        # Initialisation of graph for y-mode
        self.yGraph = pylab.Figure(figsize=(6, 4), dpi=100)
        self.yCanvas = FigureCanvasTkAgg(self.yGraph, master=self.top)
        self.yFig = self.yGraph.add_subplot(111)
        self.yFig.axes.set_xticks([])
        self.yFig.axes.set_yticks([])
        self.yFig.axis([0, 382, 255, 0])
        self.yFig.axes.set_title("Y-mode")
        self.yCanvas.show()
        self.yCanvas.get_tk_widget().grid(column=0, row=0, columnspan=3)

        # Initialisation of graph for x-mode
        self.xGraph = pylab.Figure(figsize=(6, 4), dpi=100)
        self.xCanvas = FigureCanvasTkAgg(self.xGraph, master=self.top)
        self.xFig = self.xGraph.add_subplot(111)
        self.xFig.axes.set_xticks([])
        self.xFig.axes.set_yticks([])
        self.xFig.axis([0, 382, 255, 0])
        self.xFig.axes.set_title("X-mode")
        self.xCanvas.show()
        self.xCanvas.get_tk_widget().grid(column=0, row=1, columnspan=3)

        # Add some bottons
        self.rotIncButton = Tk.Button(self.top, text="+5 deg", command=lambda: self.rotate(5.0))
        self.rotIncButton.grid(column=0, row=2, pady=10)
        self.rotDecButton = Tk.Button(self.top, text="-5 deg", command=lambda: self.rotate(-5.0))
        self.rotDecButton.grid(column=1, row=2, pady=10)
        self.okButton = Tk.Button(self.top, text="Ok", command=self.top.destroy)
        self.okButton.grid(column=2, row=2, pady=10)

    def rotate(self, angle):
        self.rotation += angle
        print self.rotation
        # Show catalogues
        self.show_cat()
        

        
    def show_cat(self):
        """ Show catalogue as an image to check if some polar pairs are too close"""
        xSize = 382
        ySize = 255
        xCoords = self.window.seCat.get_all_values("X_IMAGE")
        yCoords = self.window.seCat.get_all_values("Y_IMAGE")

        gridX, gridY = np.meshgrid(range(xSize), range(ySize))

        # Create images for x and y mode
        dataY = np.zeros((ySize, xSize))
        dataX = np.zeros((ySize, xSize))
        for obj in self.window.seCat:
            if (obj["FLUX_AUTO"] <= 0) or (obj["FWHM_IMAGE"] <= 0):
                continue
            xObj = obj["X_IMAGE"]
            yObj = obj["Y_IMAGE"]
            xPairY = xObj + 17.7
            yPairY = yObj - 0.7
            xPairX = xObj + 12.5
            yPairX = yObj - 12.5
            sqDistsObj = (gridX-xObj)**2.0 + (gridY-yObj)**2.0
            sqDistsPairY = (gridX-xPairY)**2.0 + (gridY-yPairY)**2.0
            sqDistsPairX = (gridX-xPairX)**2.0 + (gridY-yPairX)**2.0
            objImag = obj["FLUX_AUTO"] * np.exp(-sqDistsObj/(2*obj["FWHM_IMAGE"])) 
            dataY += objImag
            dataY += obj["FLUX_AUTO"] * np.exp(-sqDistsPairY/(2*obj["FWHM_IMAGE"]))
            dataX += objImag
            dataX += obj["FLUX_AUTO"] * np.exp(-sqDistsPairX/(2*obj["FWHM_IMAGE"]))

        yMean = np.mean(dataY)
        yStd = np.std(dataY)
        self.yFig.imshow(dataY, interpolation='nearest', cmap='gray', vmin=yMean, vmax=yMean+2*yStd)
        self.yCanvas.show()

        xMean = np.mean(dataX)
        xStd = np.std(dataX)
        self.xFig.imshow(dataX, interpolation='nearest', cmap='gray', vmin=xMean, vmax=xMean+2*xStd)
        self.xCanvas.show()
