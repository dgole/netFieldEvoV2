#!/usr/bin/python
from __future__ import unicode_literals
import matplotlib as m
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import math
from scipy import fftpack as fft
import sys
from matplotlib.backends.backend_pdf import PdfPages
import resource
import time
import numpy.random
import matplotlib.animation as animation
from matplotlib.pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import functionLib as lib

#m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True
font = {'family' : 'normal', 'weight' : 'bold', 'size'   : 16}
fontLabels = {'family' : 'normal', 'weight' : 'bold', 'size'   : 20}
m.rc('font', **font)

m.rcParams['xtick.major.size']  = 8
m.rcParams['xtick.major.width'] = 1.5
m.rcParams['xtick.minor.size']  = 4
m.rcParams['xtick.minor.width'] = 1

m.rcParams['ytick.major.size']  = m.rcParams['xtick.major.size']
m.rcParams['ytick.major.width'] = m.rcParams['xtick.major.width']
m.rcParams['ytick.minor.size']  = m.rcParams['xtick.minor.size']
m.rcParams['ytick.minor.width'] = m.rcParams['xtick.minor.width']

# import data as an object
# if the object "does it", make it a method, ie getting data
# if it acts on the object make it a function, ie plotting data
# functions add plots to a figure but don't export

class Data:
	def __init__(self, path, savePath=None, nRead=None):
		print("initializing data structure from " + str(path))
		self.path = path
		paramFileName = self.path+"params.txt"
		read = np.asarray(np.genfromtxt(paramFileName, dtype=str))
		paramList = []
		for i in range(read.shape[0]):
			paramList.append(read[i])
		print(read)
		print(paramList)
		self.inp  = lib.InputFile(int(paramList[0]), float(paramList[1]), float(paramList[2]), float(paramList[3]), float(paramList[4]), float(paramList[5]))
		sgrid     = np.load(self.path+"sgrid.npy")
		dgrid     = np.load(self.path+"dgrid.npy")
		state     = np.load(self.path+"state.npy")
		time      = np.load(self.path+"time.npy" )
		if nRead is not None:
			sgrid   = sgrid[0:nRead]
			dgrid   = dgrid[0:nRead]
			state   = state[0:nRead]
			time    = time[0:nRead]
		self.data = []
		for col in range(dgrid.shape[0]): self.data.append(dgrid[col])
		for col in range(state.shape[0]): self.data.append(state[col])
		if savePath is None: self.savePath = self.path
		if not os.path.exists(self.savePath): os.makedirs(self.savePath)
		if not os.path.exists(self.savePath+"/anim/"): os.makedirs(self.savePath+"/anim/")
		self.r       = sgrid[0]
		self.dr      = sgrid[1]
		self.Omega   = sgrid[2]
		self.tOrbit  = (2.0*3.14159)/self.Omega
		self.vKep    = self.Omega*self.r
		self.t       = time
		self.dt      = self.t-np.roll(self.t,1); self.dt[0]=self.dt[1]
		self.tmax    = self.t.max()
		self.nr      = self.r.shape[0]
		self.nt      = self.t.shape[0]
		self.rmax    = self.r.max()
		self.rmin    = self.r.min()
		#self.pdfName = PdfPages(self.savePath + "/plots.pdf")
		self.header  = [r"$\alpha$",
										r"$h/r$",
                    r"$T_c^4$",
                    r"$T_disk^4$",
                    r"$c_s$",
										r"$\rho$",
                    r"$\kappa_R$",
                    r"$\nu$",
                    r"$\tau$",
                    r"$\dot{M}$",
                    r"$v_{adv}$",
                    r"$v_{diff}$",
                    r"$B_z$",
                    r"$B_{rs}$",
                    r"$\psi$",
                    r"$\Sigma$",
                    r"$\beta$"]
		self.data.append(self.data[13]/self.data[12])
		self.header.append(r"$B_{rs}/B_z$")
		self.data.append(self.data[10]/self.data[11])
		self.header.append(r"$Pr_{eff}$")
		self.data.append(-self.r/self.data[10])
		self.header.append(r"$t_{adv}$")
		self.data.append(self.r/self.data[11])
		self.header.append(r"$v_{adv}$")
		self.data.append((3.0/4.0)*self.data[9]*np.square(self.Omega)*self.r*self.dr)
		self.header.append("dFlux")
		self.lum     = np.sum(self.data[21], axis=1)/np.sum(self.data[21][self.gettindex(self.inp.tWait)])
		self.lum[0]  = None
		self.logLam  = np.linspace(-10,10,num=1000)
		self.lam     = np.power(10.0, self.logLam)

	def getrindex(self, r1):
		return (np.abs(self.r-r1)).argmin()

	def gettindex(self, t1):
		return (np.abs(self.t-t1)).argmin()




class MultiPannel:
	def __init__(self, nRows=4, nCols=2, figSize=(11,13), widthRatios=[1,0.02], tmin=0, tmax=None, doExample=None):
		if tmax==None: self.tmax=doExample.tmax
		else         : self.tmax=tmax
		self.tmin=tmin
		self.nRows = nRows; self.nCols=nCols; self.figSize=figSize; 
		self.fig   = plt.figure(figsize=self.figSize)
		self.gs    = gridspec.GridSpec(nrows=nRows, ncols=2, height_ratios=np.ones(nRows), width_ratios=widthRatios)
		self.ax    = []
		self.cax   = []
		for n in range(self.nRows):
			if n==0:	self.ax.append(self.fig.add_subplot(self.gs[n, 0]))
			else   :  self.ax.append(self.fig.add_subplot(self.gs[n, 0], sharex=self.ax[0]))
			self.cax.append(self.fig.add_subplot(self.gs[n, 1]))

	def addImshowBz(self, axNum, do, tmin=0, tmax=None):
		tiMin = do.gettindex(self.tmin)
		tiMax = do.gettindex(self.tmax)
		if tmax==None: tmax=do.tmax
		extent = [self.tmin, self.tmax, np.log10(do.rmin), np.log10(do.rmax)]
		aspect = 0.3*(self.tmax-self.tmin)/(np.log10(do.rmax)-np.log10(do.rmin))
		im     = self.ax[axNum].imshow(
						    	    		    		 np.transpose(np.fliplr(do.data[12][tiMin:tiMax])),
								    		           extent=extent,
								  	  	           aspect=aspect,
										               cmap=plt.get_cmap('coolwarm'),
										               norm=colors.SymLogNorm(linthresh=0.01, linscale=1.0, vmin=-1.0, vmax=1.0)
									                 )
		self.ax[axNum].set_ylabel('log(r) (AU)')
		self.ax[axNum].set_title(do.header[12])
		self.ax[axNum].axhline(np.log10(do.inp.rIn), color='k')
		self.ax[axNum].axhline(np.log10(do.inp.rDz1), color='r')
		self.fig.colorbar(im, cax=self.cax[axNum], orientation='vertical')

	def addImshowOther(self, axNum, do, col):
		tiMin = do.gettindex(self.tmin)
		tiMax = do.gettindex(self.tmax)
		extent = [self.tmin, self.tmax, np.log10(do.rmin), np.log10(do.rmax)]
		aspect = 0.3*(self.tmax-self.tmin)/(np.log10(do.rmax)-np.log10(do.rmin))
		im     = self.ax[axNum].imshow(
						    	    		    		 np.transpose(np.fliplr(do.data[col][tiMin:tiMax])),
								    		           extent=extent,
								  	  	           aspect=aspect,
										               cmap=plt.get_cmap('viridis'),
										               norm=colors.LogNorm()
									                 )
		self.ax[axNum].set_ylabel('log(r) (AU)')
		self.ax[axNum].set_title(do.header[col])
		self.ax[axNum].axhline(np.log10(do.inp.rIn), color='k')
		self.ax[axNum].axhline(np.log10(do.inp.rDz1), color='r')
		self.fig.colorbar(im, cax=self.cax[axNum], orientation='vertical')

	def addLum(self, axNum, do):
		tiMin = do.gettindex(self.tmin)
		tiMax = do.gettindex(self.tmax)
		self.ax[axNum].semilogy(do.t[tiMin:tiMax], do.lum[tiMin:tiMax])
		self.ax[axNum].set_xlim(self.tmin, self.tmax)
		self.ax[axNum].set_xlabel('t (years)')
		self.ax[axNum].set_ylabel('Luminosity')







def timeScalesPlot(do, figNum=0):
	plt.figure(figNum)
	n1=1
	plt.loglog(do.r, do.tOrbit, label="orbital", color='k')
	plt.loglog(do.r, do.data[19][n1], label="advection", color='b')
	plt.loglog(do.r, do.data[20][n1], label="diffusion", color='r')
	#plt.loglog(do.r, do.data[20][n2], label="advection, active", color='r', linestyle='--')
	#plt.loglog(do.r, do.data[20][n2], label="diffusion, active", color='r')
	plt.xlim(do.inp.rIn, do.inp.rOut)
	plt.ylim(1.e-2, 1.e6)
	plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
	plt.ylabel("timescale (years)", fontdict=fontLabels)
	plt.xlabel("r (AU)"           , fontdict=fontLabels)

def profile(do, col, n, figNum=0):
	plt.figure(figNum)
	plt.loglog(do.r, do.data[col][n,:])
	plt.xlabel('r (AU)')
	plt.ylabel(do.header[col])


















