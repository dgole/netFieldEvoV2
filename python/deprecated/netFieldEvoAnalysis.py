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

tWait = 5.0

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
		#ntRightNow = self.data[0].shape[0]	
		#print("I want the time to be " + str(len(time)))
		#print("but I have to make it " + str(ntRightNow))
		#ntRightNow = ntRightNow - 10
		#print("but to be even safer I'm making it " + str(ntRightNow))
		self.r = sgrid[0]
		self.dr = sgrid[1]
		self.Omega = sgrid[2]
		self.tOrbit = (2.0*3.14159)/self.Omega
		self.vKep=self.Omega*self.r
		self.t = time
		self.dt = self.t-np.roll(self.t,1); self.dt[0]=self.dt[1]
		self.header = [r"$\alpha$", r"$h/r$", r"$T_c^4$", r"$T_disk^4$", r"$c_s$", r"$\rho$", r"$\kappa_R$", r"$\nu$", r"$\tau$", r"$\dot{M}$", r"$v_{adv}$", r"$v_{diff}$", r"$B_z$", r"$B_{rs}$", r"$\psi$", r"$\Sigma$", r"$\beta$"]  		
		self.pdfName = PdfPages(self.savePath + "/plots.pdf")
		self.tmax = self.t.max()
		self.nr = self.r.shape[0]
		self.nt = self.t.shape[0]
		self.data[12]*=-1.0
		self.data.append(self.data[13]/self.data[12]); self.header.append(r"$B_{rs}/B_z$")	
		self.data.append(self.data[10]/self.data[11]); self.header.append(r"$Pr_{eff}$")
		self.data.append(-self.r/self.data[10]); self.header.append(r"$t_{adv}$")		
		self.data.append(self.r/self.data[11]); self.header.append(r"$v_{adv}$")
		self.data.append((3.0/4.0)*self.data[9]*np.square(self.Omega)*self.r*self.dr); self.header.append("dFlux")
		self.lum = np.sum(self.data[21], axis=1)/np.sum(self.data[21][self.gettindex(tWait)]);
		self.logLam = np.linspace(-10,10,num=1000)
		self.lam = np.power(10.0, self.logLam)
		#self.r = sgrid[0]/100.0
		#self.dr = sgrid[1]/100.0
		self.rmax = self.r.max()
		self.rmin = self.r.min()
		#normArray = np.exp(self.t/1000.0)
		#self.lum = np.sum(self.data[21], axis=1)/np.sum(self.data[21][self.gettindex(10)])/normArray; 		
			
	def getrindex(self, r1):
		return (np.abs(self.r-r1)).argmin()

	def gettindex(self, t1):
		return (np.abs(self.t-t1)).argmin()

	def getLogSpacedTimeList(self, nProfiles):
		dtLog = np.log10(np.round(self.tmax+1.0))/(nProfiles-1.0)
		logtArray = np.asarray([n*dtLog for n in range(nProfiles)])
		tArray = np.power(10,logtArray)-1.0 
		return tArray

	def getLinSpacedTimeList(self, nProfiles):
		dt = self.tmax/(nProfiles)
		tArray = np.asarray([n*dt for n in range(nProfiles+1)])
		return tArray

	def calcSingleSpectrum(self, temp):
		spectrum = np.power(self.lam, -5)*np.power(np.exp(1.0/(self.lam*temp))-1.0, -1)
		return spectrum

	def calcSingleSpectrum2(self, time, r):
		temp = np.power(self.data[3][self.gettindex(time), self.getrindex(r)], 0.25)
		spectrum = r*self.dr[self.getrindex(r)]*np.power(self.lam, -5)*np.power(np.exp(1.0/(self.lam*temp))-1.0, -1)
		return spectrum

	#def multiSpectraPlot(self, time)


	def profile(self, col, time, vLineCoords=[], hLineCoords=[], logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None):
		n=self.gettindex(time)
		if savePath is None:
			savePath=self.savePath
		#plotData = self.data[col][n,9:179]; title = self.header[col];
		#rPlot = self.r[9:179]
		plotData = self.data[col][n,:]; title = self.header[col];
		rPlot = self.r[:]
		if logOption==0 and save is not None:
			plt.semilogx(rPlot, plotData);	plt.ylabel(self.header[col]);
		if logOption==1 and save is not None:
			plt.loglog(rPlot, np.absolute(plotData));	plt.ylabel(self.header[col]);
		if logOption==0 and save is None:
			plt.semilogx(rPlot, plotData, label=legendLabel);	plt.ylabel(self.header[col]);
		if logOption==1 and save is None:
			plt.loglog(rPlot, np.absolute(plotData), label=legendLabel);	plt.ylabel(self.header[col]);
		if ymin is None and ymax is None:
			plt.ylim(np.abs(plotData).min()*0.8, np.abs(plotData).max()*1.1)
		if ymin is not None and ymax is None:
			plt.ylim(ymin, plotData.max()*1.1)
		if ymin is None and ymax is not None:
			plt.ylim(plotData.min()*0.8, ymax)
		if ymin is not None and ymax is not None:
			plt.ylim(ymin, ymax)
		plt.xlabel("r (AU)");
		for xc in vLineCoords: plt.axvline(x=xc, color='k', linestyle='--')
		for yc in hLineCoords: plt.axhline(y=yc, color='k', linestyle='--')
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "profile_" + str(col) + "_" + str(n) + ".png", bbox_inches='tight')
			print("saved profile plot for column " + str(col) + " to png")
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print("saved profile plot for column " + str(col) + " to pdf")
		if save is None:
			plt.legend(bbox_to_anchor=(1.02, 1.0), loc=2, borderaxespad=0.)
			#plt.legend()
		if save is not None:
			plt.clf()

	def makeAnimFrames(self, cutFactor=1):
		for n in range(0, int(self.nt), cutFactor):
			print('saving anim frame for n = ' + str(n))
			fig = plt.figure()
			ax = []
			ax.append(subplot2grid((3, 2), (0, 0)))
			ax.append(subplot2grid((3, 2), (0, 1)))
			ax.append(subplot2grid((3, 2), (1, 0)))
			ax.append(subplot2grid((3, 2), (1, 1)))
			ax.append(subplot2grid((3, 2), (2, 0), colspan=2))

			rPlot = self.r[:]

			# Bz
			axNum = 0
			col = 12
			ax[axNum].plot(self.r, self.data[col][n])
			ax[axNum].set_xscale('log')
			ax[axNum].set_yscale('symlog', linthreshy=1.e-3)
			ax[axNum].set_ylim(-1.e0, 1.e0)
			ax[axNum].set_ylabel(self.header[col])
			ax[axNum].set_title(str(self.t[n]))

			# mdot
			axNum = 1
			col = 9
			ax[axNum].plot(self.r, self.data[col][n])
			ax[axNum].set_xscale('log')
			ax[axNum].set_yscale('log')
			ax[axNum].set_ylim(1.e-5, 1.e0)
			ax[axNum].set_ylabel(self.header[col])

			# alpha
			axNum = 2
			col = 0
			ax[axNum].plot(self.r, self.data[col][n])
			ax[axNum].set_xscale('log')
			ax[axNum].set_yscale('log')
			ax[axNum].set_ylim(1.e-4, 1.e-1)
			ax[axNum].set_ylabel(self.header[col])

			# Sigma
			axNum = 3
			col = 15
			ax[axNum].plot(self.r, self.data[col][n])
			ax[axNum].set_xscale('log')
			ax[axNum].set_yscale('log')
			ax[axNum].set_ylim(1.e-2, 1.e2)
			ax[axNum].set_ylabel(self.header[col])

			# luminosity
			axNum = 4
			ax[axNum].plot(self.t[0:n], self.lum[0:n])
			ax[axNum].set_xscale('linear')
			ax[axNum].set_yscale('log')
			ax[axNum].set_xlim(0,self.tmax)
			ax[axNum].set_ylim(1.e-1, 1.e2)
			ax[axNum].set_xlabel('Time (years)')
			ax[axNum].set_ylabel('Luminosity')


			for axNum in range(0,4): ax[axNum].axvline(x=self.r[78], color='k', linestyle='--')
			for axNum in [2,3]: ax[axNum].set_xlabel("r (AU)")

			plt.tight_layout()
			plt.savefig(self.savePath + "/anim/anim_" + str(n) + ".png", bbox_inches='tight')
			plt.close('all')




	def icPlots(self):
		plt.close('all')
		fig = plt.figure()
		ax = []
		ax.append(subplot2grid((2, 2), (0, 0)))
		ax.append(subplot2grid((2, 2), (0, 1)))
		ax.append(subplot2grid((2, 2), (1, 0)))
		ax.append(subplot2grid((2, 2), (1, 1)))

		axNum = 0
		col = 0
		title = self.header[col]		
		ax[axNum].loglog(self.r, self.data[col][0])
		ax[axNum].set_ylabel(title)
		ax[axNum].set_xlim([self.rmin, self.rmax])
		ax[axNum].axes.xaxis.set_ticklabels([])

		axNum = 1
		col = 15
		title = self.header[col]		
		ax[axNum].loglog(self.r, self.data[col][0])
		ax[axNum].set_ylabel(title)
		ax[axNum].set_xlim([self.rmin, self.rmax])
		ax[axNum].axes.xaxis.set_ticklabels([])

		axNum = 2
		col = 1
		title = self.header[col]		
		ax[axNum].loglog(self.r, self.data[col][0])
		ax[axNum].set_ylabel(title)
		ax[axNum].set_xlim([self.rmin, self.rmax])
		ax[axNum].set_ylim([1.e-3, 1.e0])
		ax[axNum].set_xlabel("r (AU)")

		axNum = 3
		col = 2
		title = r"$T_c$"		
		ax[axNum].loglog(self.r, np.power(self.data[col][0],0.25))
		ax[axNum].set_ylabel(title)
		ax[axNum].set_xlim([self.rmin, self.rmax])
		ax[axNum].set_xlabel("r (AU)")


		plt.tight_layout()

		plt.savefig(self.savePath + "IC.png", bbox_inches='tight', dpi=300); plt.clf()
		print("saved IC plots to png")





	def multiProfile(self, col, nProfiles, spacing="log", vLineCoords1=[], hLineCoords1=[], logOption=0, save=None, savePath=None, ymin=None, ymax=None):
		if savePath is None:
			savePath=self.savePath
		if spacing=="lin":
			tList = self.getLinSpacedTimeList(nProfiles)
		if spacing=="log":
			tList = self.getLogSpacedTimeList(nProfiles)
		if self.gettindex(tList[0])==self.gettindex(tList[1]):
			tList[1] = tList[0]+self.dt[0]
		if self.gettindex(tList[1]) == self.gettindex(tList[2]):
			tList[2] = tList[1]+self.dt[0]
		tList[nProfiles-1] = self.tmax
		for n in range(len(tList)-1): 
			self.profile(col, tList[n], logOption=logOption, legendLabel="t = " + str(np.round(self.t[self.gettindex(tList[n])],1)), ymin=ymin, ymax=ymax, vLineCoords=vLineCoords1, hLineCoords=hLineCoords1)
		n = nProfiles-1
		self.profile(col, tList[n], logOption=logOption, legendLabel="t = " + str(np.round(self.t[self.gettindex(tList[n])],-2)), ymin=ymin, ymax=ymax, vLineCoords=vLineCoords1, hLineCoords=hLineCoords1)
		if save=="png":
			plt.savefig(savePath + "multiProfile_" + str(col) + ".png", bbox_inches='tight', dpi=300); plt.clf()
			print("saved multi profile plot for column " + str(col) + " to png")
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight'); plt.clf()
			print("saved multi profile plot for column " + str(col) + " to pdf")



	def makeMultiAnim(self, timeCutFactor=10, lengthInSeconds=20, savePath=None, show=None, save=None):
		# Setup figure and subplots
		f0 = figure(num = 0, figsize = (16, 12))
	
		colList = [12, 0, 16, 15, 1, 9]		
		
		ax = []
		ax.append(subplot2grid((4, 2), (0, 0)))
		ax.append(subplot2grid((4, 2), (0, 1)))
		ax.append(subplot2grid((4, 2), (1, 0)))
		ax.append(subplot2grid((4, 2), (1, 1)))
		ax.append(subplot2grid((4, 2), (2, 0)))
		ax.append(subplot2grid((4, 2), (2, 1)))
		ax.append(subplot2grid((4, 2), (3, 0), colspan=2))
		
		i=0; ax[i].set_ylim(1.1*np.min(self.data[colList[i]]), 1.1*np.max(self.data[colList[i]]))
		i=1; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), 1.5*np.max(np.abs(self.data[colList[i]])))
		i=2; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), min(1.5*np.max(np.abs(self.data[colList[i]])),1.e18))
		i=3; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), 1.5*np.max(np.abs(self.data[colList[i]])))
		i=4; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), 1.5*np.max(np.abs(self.data[colList[i]])))
		i=5; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), 1.5*np.max(np.abs(self.data[colList[i]])))
		i=6; ax[i].set_ylim(0.7*np.min(self.lum), 1.5*np.max(self.lum))
		
		for axObject in ax[0:6] : axObject.set_xlim(self.rmin, self.rmax)
		ax[6].set_xlim(0.0, self.tmax)

		i=0; ax[i].set_ylabel(self.header[colList[i]])
		i=1; ax[i].set_ylabel(self.header[colList[i]])
		i=2; ax[i].set_ylabel(self.header[colList[i]])
		i=3; ax[i].set_ylabel(self.header[colList[i]])
		i=4; ax[i].set_ylabel(self.header[colList[i]])
		i=5; ax[i].set_ylabel(self.header[colList[i]])
		i=6; ax[i].set_ylabel(r"$L_{BB}/L_{BB,0}$")
		
		ax[6].set_xlabel("t")

		# Data Placeholders
		yp1=zeros(0)
		t=zeros(0)

		# set plots
		p0, = ax[0].semilogx(t,yp1)
		p1, = ax[1].loglog(t,yp1)
		p2, = ax[2].loglog(t,yp1)
		p3, = ax[3].loglog(t,yp1)
		p4, = ax[4].loglog(t,yp1)
		p5, = ax[5].loglog(t,yp1)
		p6, = ax[6].semilogy(t,yp1)

		def updateData(n):
			n1 = n*timeCutFactor
			i=0; p0.set_data(self.r, self.data[colList[i]][n1])
			i=1; p1.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			i=2; p2.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			i=3; p3.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			i=4; p4.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			i=5; p5.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			p6.set_data(self.t[0:n1], np.abs(self.lum[0:n1]))
			f0.suptitle("t = " + str(np.round(self.t[n1],1)), fontsize=16)
			return p0, p1, p2, p3, p4, p5, p6

		# interval: draw new frame every 'interval' ms
		# frames: number of frames to draw
		nFrames = int(self.nt/timeCutFactor)
		framesPerSecond = nFrames/lengthInSeconds
		simulation = animation.FuncAnimation(f0, updateData, blit=False, frames=nFrames, interval=10, repeat=False)
		# Uncomment the next line if you want to save the animation
		#simulation.save(filename='sim.mp4',fps=30,dpi=600)
		if show is not None:
			plt.tight_layout()
			plt.show()
		if save is not None:
			Writer = animation.writers['imagemagick']
			#Writer = animation.writers['ffmpeg']
			writer = Writer(fps=framesPerSecond, metadata=dict(artist='Me'), bitrate=1800)
			simulation.save(self.savePath+'anim1.mp4', writer=writer)



	def timeEvo(self, col, r, logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None, tStart=None, tEnd=None):
		i=self.getrindex(r) 
		if savePath is None:
			savePath=self.savePath
		if tStart is None: 
			plotData = self.data[col][...,i];
			tPlot = self.t
		else:
			tiStart = self.gettindex(tStart)
			tiEnd = self.gettindex(tEnd)
			plotData = self.data[col][tiStart:tiEnd,i]
			tPlot = self.t[0:self.gettindex(tEnd)-self.gettindex(tStart)]	
		title = self.header[col];
		if logOption==0 and save is not None:
			plt.plot(tPlot, plotData);	plt.ylabel(self.header[col]);
		if logOption==1 and save is not None:
			plt.semilogy(tPlot, np.absolute(plotData));	plt.ylabel(self.header[col]);
		if logOption==0 and save is None:
			plt.plot(tPlot, plotData, label=legendLabel);	plt.ylabel(self.header[col]);
		if logOption==1 and save is None:
			plt.semilogy(tPlot, np.absolute(plotData), label=legendLabel);	plt.ylabel(self.header[col]);
		if ymin is not None and ymax is None:
			plt.ylim(ymin, plotData.max()*1.1)
		if ymin is None and ymax is not None:
			plt.ylim(plotData.min()*0.8, ymax)
		if ymin is not None and ymax is not None:
			plt.ylim(ymin, ymax)
		plt.xlabel("t (code units)");
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "timeEvo_" + str(col) + ".png", bbox_inches='tight')
			print("saved time evo plot for column " + str(col) + " to png")
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print("saved time evo plot for column " + str(col) + " to pdf")
		if save is None:
			plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
			#plt.legend()
		if save is not None:
			plt.clf()


	def stPlot(self, col, cmapType="coolwarm", logOption=0, save=None, savePath=None, clim1=None, clim2=None, hLineCoords=[0, 1, 2], vLineCoords=None, xLabel=1, tStart=None, tEnd=None ):
		print(self.path + ": making ST plot for column " + str(col))
		if savePath is None:
			savePath=self.path
		if tStart is None: 
			plotData = self.data[col];
			tPlotMax = self.tmax
		else:
			tiStart = self.gettindex(tStart)
			tiEnd = self.gettindex(tEnd)
			plotData = self.data[col][tiStart:tiEnd]
			tPlotMax = tEnd-tStart	
		title = self.header[col];
		if logOption==1:
			plt.imshow(np.transpose(np.fliplr(plotData)), extent=[0,tPlotMax,np.log10(self.rmin),np.log10(self.rmax)], aspect=(80*tPlotMax/(self.rmax-self.rmin)), cmap=plt.get_cmap(cmapType), norm=LogNorm(vmin=clim1, vmax=clim2))
		if logOption==0:
			plt.imshow(np.transpose(np.fliplr(plotData)), extent=[0,tPlotMax,np.log10(self.rmin),np.log10(self.rmax)], aspect=(80*tPlotMax/(self.rmax-self.rmin)), cmap=plt.get_cmap(cmapType))
		if logOption==2:
			plotData=plotData/np.amax(np.abs(plotData))
			plt.imshow(np.transpose(np.fliplr(plotData)), extent=[0,tPlotMax,np.log10(self.rmin),np.log10(self.rmax)], aspect=(80*tPlotMax/(self.rmax-self.rmin)), cmap=cmapType, norm=colors.SymLogNorm(linthresh=0.0001, linscale=0.1, vmin=-1.0, vmax=1.0))
		#plt.yscale('log')
		plt.title(title)
		if xLabel == 1: plt.xlabel("Time (code units)")	
		plt.ylabel("log(r) (code units)")
		plt.colorbar(shrink=0.5) #ticks=[-10^0, -10^-1, -10^-2, -10^-3, 10^-3, 10^-2, 10^-1, 10^0])
		if (clim1 is not None and clim2 is not None):
			plt.clim(clim1,clim2)
		if (clim1 is not None and clim2 is None):
			plt.clim(clim1,np.amax(plotData))
		if (clim1 is None and clim2 is not None):
			plt.clim(np.aminx(plotData),clim2)
		if vLineCoords is not None:
			for xc in vLineCoords: plt.axvline(x=xc, color='g', linestyle=':')
		for yc in hLineCoords: plt.axhline(y=yc, color='k', linestyle='--')
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "ST_" + str(col) + ".png", bbox_inches='tight')
			print("saved ST plot for column " + str(col) + " to png")
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print("saved ST plot for column " + str(col) + " to pdf")
		plt.clf()


	def multiStPlot1(self, tFlip, savePath=None, vLinesOn=True, vLineCoords1=[50], ti1=None, ti2=None, nameTag="_"):
		print(self.path + ": making multi pannel ST plot")
		if savePath is None:
			savePath=self.path

		f, axArray = plt.subplots(4, sharex=True)

		lumRangeLogY = np.log10(np.amax(self.lum)/np.amin(self.lum))
		rangeX    = self.tmax
		
		if ti1 is None: ti1 = 0
		if ti2 is None: ti2 = self.nt-1
		
		thisNt = ti2-ti1

		#aspect1 = 0.87 *(thisNt)/self.nr
		aspect1 = (15.3*2.0) * (self.t[ti2]-self.t[ti1])/self.nr
		#aspect1 = 2.28 * (self.tmax / 30.0) * ((ti2-ti1)/ti2)

		col = 12
		axNum = 0
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2] 
		plotData=plotData/np.amax(np.abs(plotData))
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='coolwarm', 
							norm=colors.SymLogNorm(linthresh=0.001, linscale=0.5, vmin=-0.5, vmax=0.5)
							)
		#, vmin=-1.0, vmax=1.0
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel(r'$log(\frac{r}{AU})$')
		axArray[axNum].set_yticks([-1,0])		

		col = 0
		axNum = 1
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2]
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis', 
							norm=LogNorm()
							)
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel(r'$log(\frac{r}{AU})$')
		axArray[axNum].set_yticks([-1,0])
		
		col = 15
		axNum = 2
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2] 
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis',
							norm=LogNorm(), vmin=1.5*np.amin(np.absolute(plotData[:,10:])), vmax=0.3*np.amax(np.absolute(plotData[:,10:]))							
							)
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel(r'$log(\frac{r}{AU})$')
		axArray[axNum].set_yticks([-1,0])

		axNum = 3
		title = "L"		
		axArray[axNum].semilogy(self.t[ti1:ti2], self.lum[ti1:ti2])
		axArray[axNum].set_ylabel(title)
		axArray[axNum].set_xlim([self.t[ti1],self.t[ti2-1]])
		maxLog = np.log10(1.2*np.amax(np.abs(self.lum[ti1:ti2]))) 
		minLog = np.log10(0.8*np.amin(np.abs(self.lum[ti1:ti2]))) 
		axArray[axNum].set_ylim([np.power(10, minLog), np.power(10, maxLog)])
		axArray[axNum].set_xlabel("time (years)")
				
		tickVals = []
		for i in range(-5, 6):
			if minLog < i < maxLog:
				tickVals.append(np.power(10.0,i))
		axArray[axNum].set_yticks(tickVals)

		vLineCoords1 = np.zeros(int(1+(self.tmax-tWait)/(2.0*tFlip)))	
		i=0	
		while i < len(vLineCoords1):
			vLineCoords1[i] = tWait + 2*i*tFlip
			i = i+1 
		hLineCoords1 = np.asarray([-1])
		
		if vLinesOn is True:	
			for axNum in np.arange(len(axArray)):
				for xc in vLineCoords1: 
					axArray[axNum].axvline(x=xc,           color='b', linestyle=':')
					axArray[axNum].axvline(x=xc+tFlip*1.0, color='r', linestyle=':')
					axArray[axNum].axvline(x=xc+tFlip*0.5, color='k', linestyle=':')
					axArray[axNum].axvline(x=xc+tFlip*1.5, color='k', linestyle=':')
		#for axNum in [0,1,2]:
			#for yc in hLineCoords1: 
				#axArray[axNum].axhline(y=yc, color='k')

 
		plt.savefig(savePath + "MST1" + nameTag + ".png", bbox_inches='tight', dpi = 600)
		print("saved MST plot to png")
		plt.clf()



	def multiStPlot1cbar(self, tFlip, savePath=None, vLinesOn=True, vLineCoords1=[50], ti1=None, ti2=None, nameTag="_"):
		print(self.path + ": making multi pannel ST plot")
		if savePath is None:
			savePath=self.path

		#f, axArray = plt.subplots(4, sharex=True)

		lumRangeLogY = np.log10(np.amax(self.lum)/np.amin(self.lum))
		rangeX    = self.tmax
		
		if ti1 is None: ti1 = 0
		if ti2 is None: ti2 = self.nt-1
		
		thisNt = ti2-ti1

		#aspect1 = 0.87 *(thisNt)/self.nr
		aspect1 = 15.3 * (self.t[ti2]-self.t[ti1])/self.nr
		#aspect1 = 2.28 * (self.tmax / 30.0) * ((ti2-ti1)/ti2)

		col = 12
		axNum = 0
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2] 
		plotData=plotData/np.amax(np.abs(plotData))
		plt.imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='coolwarm', 
							norm=colors.SymLogNorm(linthresh=0.001, linscale=0.5, vmin=-1.0, vmax=1.0)
							)
		#axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		#axArray[axNum].set_ylabel("log(r) (AU)")
		#axArray[axNum].set_yticks([-1,0])	
		plt.colorbar(shrink=0.5);
		plt.savefig(savePath + "MST1cbar0" + nameTag + ".png", bbox_inches='tight', dpi = 600); plt.clf();	

		col = 0
		axNum = 1
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2]
		plt.imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis', 
							norm=LogNorm()
							)
		#axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		#axArray[axNum].set_ylabel("log(r) (AU)")
		#axArray[axNum].set_yticks([-1,0])
		plt.colorbar(shrink=0.5);
		plt.savefig(savePath + "MST1cbar1" + nameTag + ".png", bbox_inches='tight', dpi = 600); plt.clf();
		
		col = 15
		axNum = 2
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2] 
		plt.imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis',
							norm=LogNorm()
							)
		#axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		#axArray[axNum].set_ylabel("log(r) (AU)")
		#axArray[axNum].set_yticks([-1,0])
		plt.colorbar(shrink=0.5);
		plt.savefig(savePath + "MST1cbar2" + nameTag + ".png", bbox_inches='tight', dpi = 600); plt.clf();

 






	def multiStPlot3(self, tFlip, savePath=None, vLinesOn=True, vLineCoords1=[50], ti1=None, ti2=None, nameTag="_"):
		print(self.path + ": making multi pannel ST plot")
		if savePath is None:
			savePath=self.path

		f, axArray = plt.subplots(4, sharex=True)

		lumRangeLogY = np.log10(np.amax(self.lum)/np.amin(self.lum))
		rangeX    = self.tmax
		
		if ti1 is None: ti1 = 0
		if ti2 is None: ti2 = self.nt-1
		
		thisNt = ti2-ti1

		#aspect1 = 0.87 *(thisNt)/self.nr
		aspect1 = (2.0*15.3) * (self.t[ti2]-self.t[ti1])/self.nr
		#aspect1 = 2.28 * (self.tmax / 30.0) * ((ti2-ti1)/ti2)

		col = 12
		axNum = 0
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2] 
		plotData=plotData/np.amax(np.abs(plotData))
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='coolwarm', 
							norm=colors.SymLogNorm(linthresh=0.001, linscale=0.5, vmin=-1.0, vmax=1.0)
							)
		#, vmin=-1.0, vmax=1.0
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel("log(r) (AU)")
		axArray[axNum].set_yticks([-1,0])		

		col = 0
		axNum = 1
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2]
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis', 
							norm=LogNorm()
							)
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel("log(r) (AU)")
		axArray[axNum].set_yticks([-1,0])
		
		col = 9
		axNum = 2
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2] 
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis',
							norm=LogNorm()
							)
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel("log(r) (AU)")
		axArray[axNum].set_yticks([-1,0])

		axNum = 3
		title = "L"		
		axArray[axNum].semilogy(self.t[ti1:ti2], self.lum[ti1:ti2])
		axArray[axNum].set_ylabel(title)
		axArray[axNum].set_xlim([self.t[ti1],self.t[ti2-1]])
		maxLog = np.log10(1.1*np.amax(np.abs(self.lum[int((ti2+ti1)*0.0):ti2]))) 
		minLog = np.log10(0.9*np.amin(np.abs(self.lum[int((ti2+ti1)*0.0):ti2]))) 
		axArray[axNum].set_ylim([np.power(10, minLog), np.power(10, maxLog)])
		axArray[axNum].set_xlabel("time (years)")
				
		tickVals = []
		for i in range(-5, 6):
			if minLog < i < maxLog:
				tickVals.append(np.power(10.0,i))
		axArray[axNum].set_yticks(tickVals)

		vLineCoords1 = np.zeros(int(1+(self.tmax-tWait)/(2.0*tFlip)))	
		i=0	
		while i < len(vLineCoords1):
			vLineCoords1[i] = tWait + 2*i*tFlip
			i = i+1 
		hLineCoords1 = np.asarray([-1])
		
		if vLinesOn is True:	
			for axNum in np.arange(len(axArray)):
				for xc in vLineCoords1: 
					axArray[axNum].axvline(x=xc,           color='r', linestyle=':')
					axArray[axNum].axvline(x=xc+tFlip*1.0, color='b', linestyle=':')
					axArray[axNum].axvline(x=xc+tFlip*0.5, color='k', linestyle=':')
					axArray[axNum].axvline(x=xc+tFlip*1.5, color='k', linestyle=':')
		#for axNum in [0,1,2]:
			#for yc in hLineCoords1: 
				#axArray[axNum].axhline(y=yc, color='k')

 
		plt.savefig(savePath + "MST3" + nameTag + ".png", bbox_inches='tight', dpi = 600)
		print("saved MST plot to png")
		plt.clf()










	def multiStPlot2(self, tFlip, savePath=None):
		print(self.path + ": making multi pannel ST plot")
		if savePath is None:
			savePath=self.path

		f, axArray = plt.subplots(6, sharex=True)

		lumRangeLogY = np.log10(np.amax(self.lum)/np.amin(self.lum))
		rangeX    = self.tmax
		
		#self.lum[800] = self.lum[800]*10

		#ti1 = 3000
		#ti2 = 4000
		ti1 = 0
		ti2 = self.nt
		
		#aspect1 = 1.65 * (30.0 / self.tmax)
		aspect1 = 1.49 * (self.tmax / 30.0) * ((ti2-ti1)/ti2)


		col = 12
		axNum = 0
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2] 
		plotData=plotData/np.amax(np.abs(plotData))
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='coolwarm', 
							norm=colors.SymLogNorm(linthresh=0.0001, linscale=0.1, vmin=-1.0, vmax=1.0)
							)
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel("log(r)")
		axArray[axNum].set_yticks([-1,0])		

		col = 0
		axNum = 1
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2]
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis', 
							norm=LogNorm()
							)
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel("log(r)")
		axArray[axNum].set_yticks([-1,0])
		
		col = 15
		axNum = 2
		title = self.header[col]		
		plotData = self.data[col][ti1:ti2] 
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis',
							norm=LogNorm()
							)
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel("log(r)")
		axArray[axNum].set_yticks([-1,0])

		col = 2
		axNum = 3
		title = self.header[col]		
		plotData = np.power(self.data[col][ti1:ti2], 1.0) 
		axArray[axNum].imshow(np.transpose(np.fliplr(plotData)), 
							extent=[self.t[ti1],self.t[ti2-1],np.log10(self.rmin),np.log10(self.rmax)], 
							aspect=aspect1, 
							cmap='viridis',
							norm=LogNorm()
							)
		axArray[axNum].text(1.02, 0.5, title, fontsize=14,transform=axArray[axNum].transAxes)
		axArray[axNum].set_ylabel("log(r)")
		axArray[axNum].set_yticks([-1,0])

		axNum = 4
		title = "L"		
		axArray[axNum].semilogy(self.t[ti1:ti2], self.lum[ti1:ti2])
		axArray[axNum].set_ylabel(title)
		axArray[axNum].set_xlim([self.t[ti1],self.t[ti2-1]])
		maxLog = np.log10(1.1*np.amax(np.abs(self.lum[int((ti2+ti1)*0.0):ti2]))) 
		minLog = np.log10(0.9*np.amin(np.abs(self.lum[int((ti2+ti1)*0.0):ti2]))) 
		axArray[axNum].set_ylim([np.power(10, minLog), np.power(10, maxLog)])
				
		tickVals = []
		for i in range(-5, 6):
			if minLog < i < maxLog:
				tickVals.append(np.power(50.0,i))
		axArray[axNum].set_yticks(tickVals)

		axNum = 5
		title = self.header[9]		
		plotData = self.data[9][ti1:ti2,10]/self.data[9][self.gettindex(tWait),10]
		axArray[axNum].semilogy(self.t[ti1:ti2], plotData)
		axArray[axNum].set_ylabel(title)
		axArray[axNum].set_xlim([self.t[ti1],self.t[ti2-1]])
		axArray[axNum].set_xlabel("Time")
		maxLog = np.log10(1.1*np.amax(np.abs(plotData[int((ti2+ti1)*0.0):ti2]))) 
		minLog = np.log10(0.9*np.amin(np.abs(plotData[int((ti2+ti1)*0.0):ti2]))) 
		axArray[axNum].set_ylim([np.power(10, minLog), np.power(10, maxLog)])
		tickVals = []
		for i in range(-5, 6):
			if minLog < i < maxLog:
				tickVals.append(np.power(10.0,i))
		axArray[axNum].set_yticks(tickVals)

		#vLineCoords1 = np.zeros(int(1+(self.tmax-50.0)/(2.0*tFlip)))	
		#i=0	
		#while i < len(vLineCoords1):
			#vLineCoords1[i] = 50.0 + 2*i*tFlip
			#i = i+1 
		#hLineCoords1 = np.asarray([-1])

		#if self.data[12][self.gettindex(50.0), 5] > self.data[12][self.gettindex(10.5), 5]:
		#	color1 = 'b'
		#	color2 = 'r'
		#else:
		#	color1 = 'r'
		#	color2 = 'b'			
		
		#for axNum in np.arange(len(axArray)):
		#	for xc in vLineCoords1: 
		#		axArray[axNum].axvline(x=xc, color=color1, linestyle=':')
		#		axArray[axNum].axvline(x=xc+tFlip, color=color2, linestyle=':')
		#for axNum in [0,1,2]:
			#for yc in hLineCoords1: 
				#axArray[axNum].axhline(y=yc, color='k')

 
		plt.savefig(savePath + "MST2.png", bbox_inches='tight', dpi = 600)
		print("saved MST plot to png")
		plt.clf()

		


'''
	def stPlot(self, col, cmapType="coolwarm", logOption=0, save=None, savePath=None, clim1=None, clim2=None):
		print(self.path + ": making ST plot for column " + str(col))
		if savePath is None:
			savePath=self.path
		if logOption==0:
			plotData = self.data[col]; title = self.header[col];
		if logOption==1:
			plotData = np.log10(np.absolute(self.data[col])); title = "log " + self.header[col];
		plt.imshow(np.transpose(np.fliplr(plotData)), extent=[0,self.tmax,self.rmin,self.rmax], aspect=(0.5*self.tmax/(self.rmax-self.rmin)), cmap=plt.get_cmap(cmapType))
		plt.yscale('log')
		plt.title(title); plt.xlabel("Time (code units)"); plt.ylabel("r (code units)");
		plt.colorbar(shrink=0.5)
		if (clim1 is not None and clim2 is not None):
			plt.clim(clim1,clim2)
		if (clim1 is not None and clim2 is None):
			plt.clim(clim1,np.amax(plotData))
		if (clim1 is None and clim2 is not None):
			plt.clim(np.aminx(plotData),clim2)
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "ST_" + str(col) + ".png", bbox_inches='tight')
			print("saved ST plot for column " + str(col) + " to png")
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print("saved ST plot for column " + str(col) + " to pdf")
		plt.clf()
'''




'''
	def addCol(self, funcName, label, *args, **kwargs):
		print(self.path + ": adding column named " + label)
		self.data.append(funcName(self.data))
		self.header.append(label)

	def multiTimeEvo(self, col, rList, logOption=0, save=None, savePath=None, ymin=None, ymax=None):
		for r in rList: 
			self.timeEvo(col, r, logOption=logOption, legendLabel="r = " + str(np.round(r,1)), ymin=ymin, ymax=ymax)
		if save=="png":
			plt.savefig(savePath + "multiTimeEvo_" + str(col) + ".png", bbox_inches='tight')
			print "saved multi time evo plot for column " + str(col) + " to png"
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print "saved multi time evo plot for column " + str(col) + " to pdf"
		plt.clf()

	def makeAnim(self, animCol, logOption=0, savePath=None, show=None, save=None):
		if savePath is None:
			savePath=self.savePath
		fig = plt.figure()
		if logOption==1:
			ylim1=np.amin(np.abs(self.data[animCol]))*0.8
			ylim2=np.amax(np.abs(self.data[animCol]))*1.1
		if logOption==0:
			ylim1=np.amin(self.data[animCol])
			ylim2=np.amax(self.data[animCol])
			if ylim1<0:
				ylim1=ylim1*1.1
			else:
				ylim1=ylim1*0.8
			if ylim2<0:
				ylim2=ylim2*1.1
			else:
				ylim2=ylim2*0.8
		ax = plt.axes(xlim=(self.r[0], self.r[self.r.shape[0]-1]), ylim=(ylim1, ylim2))
		if logOption==1:
			line, = ax.loglog([], [], lw=2)
		if logOption==0:
			line, = ax.semilogx([], [], lw=2)
		def init():
			line.set_data([], [])
			return line, 
		if logOption==1:
			def animate(i):
				x = self.r
				y = np.abs(self.data[animCol][i,...])
				line.set_data(x, y)
				return line,
		if logOption==0:
			def animate(i):
				x = self.r
				y = self.data[animCol][i,...]
				line.set_data(x, y)
				return line,
		anim = animation.FuncAnimation(fig, animate, init_func=init, frames=self.nt, interval=20, blit=True)
		plt.xlabel("r (code units)")
		plt.ylabel(self.header[animCol])
		if show is not None:
			plt.tight_layout()
			plt.show()
		if save is not None:
			Writer = animation.writers['ffmpeg']
			writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
			anim.save(savePath+'anim.mp4', writer=writer)

	def stPlot(self, col, cmapType="viridis", logOption=0, save=None, savePath=None, clim1=None, clim2=None):
		print self.path + ": making ST plot for column " + str(col)
		if savePath is None:
			savePath=self.path
		if logOption==0:
			plotData = self.data[col]; title = self.header[col];
		if logOption==1:
			plotData = np.log10(np.absolute(self.data[col])); title = "log " + self.header[col];
		plt.imshow(np.transpose(np.fliplr(plotData)), extent=[0,self.tmax,self.rmin,self.rmax], aspect=(0.2*self.tmax/(self.rmax-self.rmin)), cmap=plt.get_cmap(cmapType))
		plt.title(title); plt.xlabel("Time (code units)"); plt.ylabel("r (code units)");
		plt.colorbar(shrink=0.5)
		if (clim1 is not None and clim2 is not None):
			plt.clim(clim1,clim2)
		if (clim1 is not None and clim2 is None):
			plt.clim(clim1,np.amax(plotData))
		if (clim1 is None and clim2 is not None):
			plt.clim(np.aminx(plotData),clim2)
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "ST_" + str(col) + ".png", bbox_inches='tight')
			print "saved ST plot for column " + str(col) + " to png"
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print "saved ST plot for column " + str(col) + " to pdf"
		plt.clf()

	def addCol(self, funcName, label, *args, **kwargs):
		print self.path + ": adding column named " + label
		self.data.append(funcName(self.data))
		self.header.append(label)


		
'''











