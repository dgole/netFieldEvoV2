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
import matplotlib.colors as colors
m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True

font = {'family' : 'normal', 'weight' : 'bold', 'size'   : 10}
m.rc('font', **font)

class Data:
	def getrindex(self, r1):
		return (np.abs(self.r-r1)).argmin()
	def gettindex(self, t1):
		return (np.abs(self.t-t1)).argmin()
	def __init__(self, path, savePath=None, nRead=None):
		print("initializing data structure from " + str(path))
		self.path= path
		sgrid = np.load(self.path+"sgrid.npy")
		dgrid = np.load(self.path+"dgrid.npy")
		state = np.load(self.path+"state.npy")
		time  = np.load(self.path+"time.npy")
		if nRead is not None:
			sgrid = sgrid[0:nRead]
			dgrid = dgrid[0:nRead]
			state = state[0:nRead]
			time  = time[0:nRead]
		self.data = []
		for col in range(dgrid.shape[0]): self.data.append(dgrid[col])
		for col in range(state.shape[0]): self.data.append(state[col])
		if savePath is None:
			self.savePath = self.path
		if not os.path.exists(self.savePath): os.makedirs(self.savePath)
		if not os.path.exists(self.savePath+"/anim/"): os.makedirs(self.savePath+"/anim/")
		self.r = sgrid[0]
		self.dr = sgrid[1]
		self.Omega = sgrid[2]
		self.tOrbit = (2.0*3.14159)/self.Omega
		self.vKep=self.Omega*self.r
		self.t = time
		self.dt = self.t-np.roll(self.t,1); self.dt[0]=self.dt[1]
		self.header = [r"$\alpha$", r"$h/r$", r"$T_c^4$", r"$T_disk^4$", r"$c_s$", r"$\rho$", r"$\kappa_R$", r"$\nu$", r"$\tau$", r"$\dot{M}$", r"$-v_{adv}$", r"$v_{diff}$", r"$B_z$", r"$B_{rs}$", r"$\psi$", r"$\Sigma$", r"$\beta$"]  		
		self.tmax = self.t.max()
		self.nr = self.r.shape[0]
		self.nt = self.t.shape[0]
		self.data[10]*=-1.0
		self.data.append(self.data[13]/self.data[12]); self.header.append(r"$B_{rs}/B_z$")	
		self.data.append(self.data[10]/self.data[11]); self.header.append(r"$Pr_{eff}$")
		self.data.append(self.r/self.data[10]); self.header.append(r"$t_{adv}$")		
		self.data.append(self.r/self.data[11]); self.header.append(r"$t_{diff}$")
		self.data.append((3.0/4.0)*self.data[9]*np.square(self.Omega)*self.r*self.dr); self.header.append("dFlux")
		self.lum = np.sum(self.data[21], axis=1)/np.sum(self.data[21][self.gettindex(10)]);
		self.logLam = np.linspace(-10,10,num=1000)
		self.lam = np.power(10.0, self.logLam)
		self.rmax = self.r.max()
		self.rmin = self.r.min()

def getRGBlist(nTot):
	rgbList = [];
	for n in range(nTot): 
		r = n/nTot
		g = 0.0
		b = 1.0-n/nTot 
		rgbList.append((r,g,b))
	return rgbList	

def getLogSpacedTimeList(do, nProfiles):
	dtLog = np.log10(np.round(do.tmax+1.0))/(nProfiles-1.0)
	logtArray = np.asarray([n*dtLog for n in range(nProfiles)])
	tArray = np.power(10,logtArray)-1.0 
	return tArray

def getLinSpacedTimeList(do, nProfiles):
	dt = do.tmax/(nProfiles)
	tArray = np.asarray([n*dt for n in range(nProfiles+1)])
	return tArray

def getProfileData(do, col, t):
	return do.r, do.data[col][do.gettindex(t)]

def profilePlot(do, col, t, figNum=0, logOption='loglog', color='k', label='time'):
	plt.figure(0)
	x, y = getProfileData(do, col, t)
	if logOption=='loglog':	
		if label=='time':
			plt.loglog(x, y, color=color, label="t="+str(np.round(t,1)))
		elif label=='header':  
			plt.loglog(x, y, color=color, label=do.header[col])
	if logOption=='semilogx':	
		if label=='time':
			plt.semilogx(x, y, color=color, label="t="+str(np.round(t,1)))
		elif label=='header':  
			plt.semilogx(x, y, color=color, label=do.header[col])
	plt.xlabel(r'$r$'); plt.ylabel(do.header[col]);

def multiProfilePlot(do, col, nProfiles, tOption='log', figNum=0, logOption='loglog'):
	if   tOption == 'log':
		tList = getLogSpacedTimeList(do, nProfiles)
	elif tOption == 'lin':
		tList = getLinSpacedTimeList(do, nProfiles)
	rgbList = getRGBlist(nProfiles)
	n=0;
	for t in tList:
		profilePlot(do, col, t, logOption=logOption, color=rgbList[n])
		n+=1;
	plt.legend(loc=(1.01, 0.0))
		








