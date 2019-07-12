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
import netFieldEvoAnalysis as reader
import resource
from matplotlib.backends.backend_pdf import PdfPages
import time
import matplotlib.animation as animation
m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True

print(sys.path)
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

# avaliable routines
# profile(self, col, time, vLineCoords=[1.0, 50.0], hLineCoords=[], logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None)
# multiProfile(self, col, nProfiles, spacing="log", vLineCoords1=[1.0,50.0], hLineCoords1=[], logOption=0, save=None, savePath=None, ymin=None, ymax=None)
# makeMultiAnim(self, timeCutFactor=10, lengthInSeconds=20, savePath=None, show=None, save=None)
# timeEvo(self, col, r, logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None)
# stPlot(self, col, cmapType="viridis", logOption=0, save=None, savePath=None, clim1=None, clim2=None)

for idNum in sys.argv[2:]:
	# make data object
	do = reader.Data("../output/run" + str(idNum) + "/")

	for i in range(do.nr):
		print(i, do.r[i], do.data[20][0,i], do.data[19][0,i])

	print(do.t)

	#do.multiStPlot(float(sys.argv[1]))

	#print(do.data[12][3800,90]/do.data[12][0,90])
	#print(do.data[12][3800,90]/do.data[12][0,90])
	#print(do.data[12][3800,130]/do.data[12][0,130])
	#print(do.data[12][3800,130]/do.data[12][0,130])

	# all time scales on the same plot
	n=0
	plt.loglog(do.r[9:179], do.tOrbit[9:179], label="Orbit")
	plt.loglog(do.r[9:179], do.data[19][n,9:179], label="Advection")
	plt.loglog(do.r[9:179], do.data[20][n,9:179], label="Diffusion")
	plt.legend(bbox_to_anchor=(0.7, 0.3), loc=2, borderaxespad=0.)
	plt.ylim(1.e-2,1.e3)
	plt.ylabel("Time-scale (years)")
	plt.xlabel("r (AU)")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight')
	plt.savefig(do.savePath + "timeScales.png", bbox_inches='tight', dpi=300)
	plt.clf()

	# space time plot
	do.stPlot(12, cmapType="coolwarm", logOption=2, save="pdf")
	do.stPlot(0, cmapType="viridis", logOption=1, save="pdf")
	do.stPlot(9, cmapType="viridis", logOption=1, save="pdf", clim1=1.e8, clim2=1.e10)

	do.stPlot(12, cmapType="coolwarm", logOption=2, save="png", xLabel=0)
	do.stPlot(0, cmapType="viridis", logOption=1, save="png", xLabel=1)
	do.stPlot(9, cmapType="viridis", logOption=1, save="png", clim1=1.e8, clim2=1.e10, xLabel=1)

	# sigma multi profile
	do.multiProfile(15, 10, logOption=1, save="pdf", spacing="log")
	# alpha multi profile
	do.multiProfile(0, 10, logOption=1, save="pdf", spacing="log")
	# h/r multi profile
	do.multiProfile(1, 10, logOption=1, ymin=1.e-2, ymax=1.e0, save="pdf", spacing="log")
	# nu multi profile
	do.multiProfile(7, 10, logOption=1, save="pdf", spacing="log")
	# mdot multi profile
	do.multiProfile(9, 10, logOption=1, save="pdf", spacing="log")
	# bz multi profile
	do.multiProfile(12, 10, logOption=0, save="pdf", spacing="log")
	do.multiProfile(12, 10, logOption=1, save="png", ymin=5.e-4, ymax=7e0)
	# beta multi profile
	do.multiProfile(16, 10, logOption=1, save="pdf")

	# brs multi profile
	#do.multiProfile(6, 10, logOption=0, save="pdf", ymin=-1.0, ymax=1.0)
	#do.multiProfile(6, 10, logOption=1, save="pdf", ymin=1.e-3, ymax=1.e3)
	# psi multi profile
	do.multiProfile(14, 10, logOption=1, spacing="log", save="pdf")
	# brs/bz multi profile
	do.multiProfile(17, 10, logOption=0, save="pdf", ymin=-2.0, ymax=2.0)




	do.profile(0, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	do.profile(1, 0.01, vLineCoords=[10.0, 100.0], logOption=1, ymin=1.e-2, ymax=1.e0, save='png')
	do.profile(9, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	do.profile(12, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	do.profile(15, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')

	do.pdfName.close()
	#tStart1 = 10; tEnd1 = 32*10;
'''
	tStart1 = None; tEnd1 = None;
	dt = float(sys.argv[1])
	vLineCoords1 = np.zeros(int(1+(do.tmax-10.0)/dt))
	vLineCoords2 = np.zeros(int(1+(do.tmax-10.0)/dt))
	i=0	
	print(int(1+(do.tmax-10.0)/dt))
	while i < len(vLineCoords1):
		vLineCoords1[i] = 10.0 + i*dt
		i = i+1 
	print(dt)
	print(vLineCoords1)
	#vLineCoords1 = None	
	do.stPlot(21, cmapType="viridis", logOption=1, save="png", xLabel=1, tStart=tStart1, tEnd=tEnd1, vLineCoords=vLineCoords1)
	do.stPlot(12, cmapType="coolwarm", logOption=2, save="png", xLabel=1, tStart=tStart1, tEnd=tEnd1, vLineCoords=vLineCoords1)
	do.stPlot(15, cmapType="viridis", logOption=1, save="png", xLabel=1, tStart=tStart1, tEnd=tEnd1, vLineCoords=vLineCoords1)
	do.stPlot(0, cmapType="viridis", logOption=1, save="png", xLabel=1, tStart=tStart1, tEnd=tEnd1, vLineCoords=vLineCoords1)
	do.stPlot(9, cmapType="viridis", logOption=1, save="png", clim1=1.e8, clim2=1.e10, xLabel=1, tStart=tStart1, tEnd=tEnd1, vLineCoords=vLineCoords1)
	do.timeEvo(9, 3.0, logOption=1, save='png', tStart=tStart1, tEnd=tEnd1)

	plt.figure(figsize=(10, 4))
	print(do.lum.shape)
	if tStart1 is not None:
		plt.semilogy(do.t[0:do.gettindex(tEnd1)-do.gettindex(tStart1)], do.lum[do.gettindex(tStart1):do.gettindex(tEnd1)])
	else:
		plt.semilogy(do.t, do.lum)
	plt.ylabel("Integrated Luminosity (code units)")
	plt.xlabel("Time (code units)")
	if vLineCoords1 is not None: 
		for xc in vLineCoords1: 
			plt.axvline(x=xc, color='g', linestyle=':')
	#plt.axhline(y=1.e0)
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight')
	plt.savefig(do.savePath + "lightCurve.png", bbox_inches='tight')
	plt.clf()
'''
	
	#print(do.lum[-1]/do.lum[-100])
	











'''
if str(sys.argv[2])=="save": 
	t0 = time.time()
	do.makeMultiAnim(save="yes", timeCutFactor=int(sys.argv[3]), lengthInSeconds=int(60))
	print(time.time() - t0)
if str(sys.argv[2])=="show": do.makeMultiAnim(show="yes")
'''


























