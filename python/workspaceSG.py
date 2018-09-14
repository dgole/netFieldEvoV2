#!/usr/bin/python
from __future__ import unicode_literals
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import math
from scipy import fftpack as fft
import sys
import netFieldEvoAnalysisSG as reader
import resource
from matplotlib.backends.backend_pdf import PdfPages
import time
import matplotlib.animation as animation
m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True

######################
#self.header = [r"$\alpha$", r"$h/r$", r"$T_c^4$", r"$T_disk^4$", r"$c_s$", r"$\rho$", r"$\kappa_R$", r"$\nu$", r"$\tau$", r"$mdot$", r"$v_{adv}$", r"$v_{diff}$", r"$\beta$", r"$B_z$", r"$B_{rs}$", r"$\psi$", r"$\Sigma$"  		
		
# 0  alpha
# 1  h/r
# 2  Tc4
# 3  Tdisk4
# 4  cs
# 5  rho
# 6  kappaR
# 7  nu
# 8  tau
# 9  mdot
# 10 vadv
# 11 vdiff
# 12 Bz
# 13 Brs
# 14 Psi
# 15 Sigma
# 16 beta
# 17 Bz/Brs
# 18 Pr eff
# 19 tadv 
# 20 tdiff
# 21 dFlux
# 22 Qinv

# avaliable routines
# profile(self, col, time, vLineCoords=[1.0, 50.0], hLineCoords=[], logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None)
# multiProfile(self, col, nProfiles, spacing="log", vLineCoords1=[1.0,50.0], hLineCoords1=[], logOption=0, save=None, savePath=None, ymin=None, ymax=None)
# makeMultiAnim(self, timeCutFactor=10, lengthInSeconds=20, savePath=None, show=None, save=None)
# timeEvo(self, col, r, logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None)
# stPlot(self, col, cmapType="viridis", logOption=0, save=None, savePath=None, clim1=None, clim2=None)

for idNum in sys.argv[1:]:
	# make data object
	do = reader.Data("../outputSG/run" + str(idNum) + "/")

	for i in range(do.nr):
		print(i, do.r[i], do.tOrbit[i], do.data[20][0,i], do.data[19][0,i])

	print(do.data[1].shape)
	print(do.data[10][100])
	print(do.data[11][100])

	do.icPlots()
	do.multiStPlot1(float(sys.argv[1]))
	do.multiStPlot2(float(sys.argv[1]))
	

	# all time scales on the same plot
	n=0
	plt.loglog(do.r, do.tOrbit, label="tOrbit", color='k', linestyle='-')
	plt.loglog(do.r, do.data[19][n], label="tAdv", color='b', linestyle='-')
	plt.loglog(do.r, do.data[20][n], label="tDiff", color='r', linestyle='-')
	n=do.nt-1
	plt.loglog(do.r, do.data[19][n], label="tAdv late", color='b', linestyle='--')
	plt.loglog(do.r, do.data[20][n], label="tDiff late", color='r', linestyle='--')
	plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
	plt.ylim(1.e-1,1.e8)
	plt.grid()
	plt.ylabel("timescale (years)")
	plt.xlabel("r (AU)")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight')
	plt.savefig(do.savePath + "timeScales.png", bbox_inches='tight')
	plt.clf()

	# space time plot
	#do.stPlot(12, cmapType="coolwarm", logOption=2, save="pdf")
	#do.stPlot(0, cmapType="viridis", logOption=1, save="pdf")
	#do.stPlot(9, cmapType="viridis", logOption=1, save="pdf", clim1=1.e8, clim2=1.e10)

	#do.stPlot(12, cmapType="coolwarm", logOption=2, save="png", xLabel=0)
	#do.stPlot(0, cmapType="viridis", logOption=1, save="png", xLabel=1)
	#do.stPlot(9, cmapType="viridis", logOption=1, save="png", clim1=1.e8, clim2=1.e10, xLabel=1)

	# sigma multi profile
	do.multiProfile(15, 10, logOption=1, ymin=np.amin(do.data[15]), ymax=np.amax(do.data[15]), save="pdf", spacing="log")
	do.multiProfile(15, 10, logOption=1, ymin=np.amin(do.data[15]), ymax=np.amax(do.data[15]), save="pdf", spacing="lin")
	# alpha multi profile
	do.multiProfile(0, 10, logOption=1, ymin=0.9*np.amin(do.data[0]), ymax=1.1*np.amax(do.data[0]), save="pdf", spacing="log")
	do.multiProfile(0, 10, logOption=1, ymin=0.9*np.amin(do.data[0]), ymax=1.1*np.amax(do.data[0]), save="pdf", spacing="lin")
	# Qinv multi profile
	do.multiProfile(22, 10, logOption=1, ymin=np.amin(do.data[22]), ymax=np.amax(do.data[22]), spacing="log", save="pdf", hLineCoords1=[1.0])
	# h/r multi profile
	do.multiProfile(1, 10, logOption=1, ymin=0.5*np.amin(do.data[1]), ymax=2*np.amax(do.data[1]), save="pdf", spacing="log")
	# nu multi profile
	do.multiProfile(7, 10, logOption=1, ymin=np.amin(do.data[7]), ymax=np.amax(do.data[7]), save="pdf", spacing="log")
	# mdot multi profile
	do.multiProfile(9, 10, logOption=1, ymin=np.amin(do.data[9]), ymax=np.amax(do.data[9]), save="pdf", spacing="log")
	# bz multi profile
	do.multiProfile(12, 10, logOption=0, ymin=np.amin(do.data[12]), ymax=np.amax(do.data[12]), save="pdf", spacing="log")
	do.multiProfile(12, 10, logOption=0, ymin=np.amin(do.data[12]), ymax=np.amax(do.data[12]), save="pdf", spacing="lin")
	do.multiProfile(12, 10, logOption=1, ymin=0.9*np.amin(do.data[12]), ymax=1.1*np.amax(do.data[12]), save="pdf", spacing="log")
	do.multiProfile(12, 10, logOption=1, ymin=0.9*np.amin(do.data[12]), ymax=1.1*np.amax(do.data[12]), save="pdf", spacing="lin")
	# beta multi profile
	do.multiProfile(16, 10, logOption=1, ymin=np.amin(do.data[16]), ymax=np.amax(do.data[16]), save="pdf")

	# brs multi profile
	do.multiProfile(13, 10, logOption=0, save="pdf", ymin=-1.0, ymax=1.0)
	do.multiProfile(13, 10, logOption=1, save="pdf", ymin=np.amin(do.data[13]), ymax=np.amax(do.data[13]))
	# psi multi profile
	do.multiProfile(14, 10, logOption=1, ymin=np.amin(do.data[14]), ymax=np.amax(do.data[14]), spacing="log", save="pdf")
	# brs/bz multi profile
	do.multiProfile(17, 10, logOption=0, save="pdf", ymin=-1.0, ymax=1.0)
	
	do.multiProfile(2, 10, logOption=1, ymin=np.amin(do.data[2]), ymax=np.amax(do.data[2]), spacing="log", save="pdf")


	
	plt.loglog(do.t, do.mtot)
	plt.xlabel("time (years)"); plt.ylabel("total mass")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight'); plt.clf();

	plt.plot(do.t, do.Bztot)
	plt.xlabel("time (years)"); plt.ylabel("total vertical magnetic flux")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight'); plt.clf();

	plt.loglog(do.t, do.Bztot)
	plt.xlabel("time (years)"); plt.ylabel("total vertical magnetic flux")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight'); plt.clf();

	plt.plot(do.t, do.Bztot2)
	plt.xlabel("time (years)"); plt.ylabel("normal domain vertical magnetic flux")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight'); plt.clf();

	plt.loglog(do.t, do.Bztot2)
	plt.xlabel("time (years)"); plt.ylabel("normal domain vertical magnetic flux")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight'); plt.clf();


	transLoc = np.zeros_like(do.t)
	for n in range(do.nt):
		transLoc[n] = do.r[np.abs(do.data[22][n,0:180]-1.0).argmin()] 
	print(transLoc)
	plt.loglog(do.t, transLoc)
	plt.xlabel("time (years)"); plt.ylabel("SG transition location (AU)")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight'); plt.clf();



	

	do.pdfName.close()



