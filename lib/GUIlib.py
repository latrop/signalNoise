#! /usr/bin/env python

from os import path
import Tkinter as Tk
import tkFileDialog, tkMessageBox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import pylab
import pyfits
import numpy as np

from lib.alignment import *
from lib.se import clean_background


class MenuBar(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.menubar = Tk.Menu(window.root)
        self.fileMenu = Tk.Menu(self.menubar, tearoff=0)
        # self.fileMenu.add_command(label="Select object", command=self.select_object)
        self.fileMenu.add_command(label="Select folder", command=self.select_folder)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
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
        dirPath = tkFileDialog.askdirectory(parent=self.window.root,
                                            title="Open data folder")
        dirPath = path.normpath(dirPath)
        self.window.setup(dirPath)


class ImagPanel(Tk.Frame):
    def __init__(self, window):
        self.window = window
        self.mainGraph = pylab.Figure(figsize=(6, 4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.mainGraph, master=self.window.root)
        self.fig = self.mainGraph.add_subplot(111)
        self.fig.axes.set_xticks([])
        self.fig.axes.set_yticks([])
        self.objPlotInstance = None
        self.standartsPlotIntance = None
        self.objPairPlotInstance = None
        self.standartsPairPlotInstance = None

        # self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.window.root)
        # self.toolbar.update()
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=0, row=0) #pack(side=Tk.TOP)

    def remove_objects_from_plot(self, upDateFigure=False):
        if self.objPlotInstance:
            self.objPlotInstance.remove()
            self.objPlotInstance = None
        if self.standartsPlotIntance:
            self.standartsPlotIntance.remove()
            self.standartsPlotIntance = None
        if self.objPairPlotInstance:
            self.objPairPlotInstance.remove()
            self.objPairPlotInstance = None
        if self.standartsPairPlotInstance:
            self.standartsPairPlotInstance.remove()
            self.standartsPairPlotInstance = None
        if upDateFigure:
            self.canvas.show()

    def plot_fits_file(self, fitsName):
        hdu = pyfits.open(fitsName)
        data = hdu[0].data.copy()
        hdu.close()
        ySize, xSize = data.shape
        meanValue = np.mean(data)
        stdValue = np.std(data)
        self.fig.imshow(data, interpolation='nearest', cmap='gray',
                       vmin=meanValue, vmax=meanValue+2*stdValue)
        self.fig.axis([0, xSize, 0, ySize])
        self.canvas.show()
        #self.canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

    def plot_objects(self, reference, polarMode=None):
        self.remove_objects_from_plot()
        self.objPlotInstance = self.fig.plot([reference.objSEParams.xImage], [reference.objSEParams.yImage], marker="o",
                                             markerfacecolor="none", markersize=15, markeredgewidth=2,
                                             markeredgecolor="r")[0]
        stx = [st['seParams'].xImage-1 for st in reference.standartsObs]
        sty = [st['seParams'].yImage-1 for st in reference.standartsObs]
        self.standartsPlotIntance =  self.fig.plot(stx, sty, marker="o", markerfacecolor="none",
                                                   linestyle="", markersize=15, markeredgewidth=2,
                                                   markeredgecolor="g")[0]

        if polarMode:
            self.objPairPlotInstance = self.fig.plot([reference.objPairSEParams.xImage],
                                                     [reference.objPairSEParams.yImage],
                                                     marker="o", markerfacecolor="none", markersize=15,
                                                     markeredgewidth=2, markeredgecolor="r")[0]
            stx = [st['seParams'].xImage-1 for st in reference.standartPairsObs]
            sty = [st['seParams'].yImage-1 for st in reference.standartPairsObs]
            self.standartsPairPlotInstance =  self.fig.plot(stx, sty, marker="o", markerfacecolor="none",
                                                            linestyle="", markersize=15, markeredgewidth=2,
                                                            markeredgecolor="g")[0]

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
        
    def update_object_info(self, objName, filtName, addString):
        s = "Object: %s" % (objName)
        if addString:
            s += "(%s)" % (addString)
        s += ", filter %s" % (filtName)
        self.objectInfoLabelValue.set(s)

    def update_message(self, key, msg):
        self.messages[key] = msg
        self.messagesString.set("\n".join(["%s: %s" % (key, self.messages[key]) for key in self.messages if self.messages[key]]))

    def show_photometry_data(self, objSn, objMag, stSn):
        string = "Mag: %1.2f\n" % (objMag)
        string += "Obj S/N: %1.0f\n" % (objSn)
        string += "\n".join(["%s S/N: %1.0f" % (st[0], st[1]) for st in stSn])
        self.photometryString.set(string)

    def show_photometry_data_polar_mode(self, objSN, objPairSN, stSnList):
        string = "Obj S/N: %1.0f/%1.0f\n" % (objSN, objPairSN)
        string += "\n".join(["%s S/N: %1.0f/%1.0f" % (st[0], st[1], st[2]) for st in stSnList])
        self.photometryString.set(string)


