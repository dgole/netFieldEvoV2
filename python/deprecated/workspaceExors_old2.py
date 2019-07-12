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
	do = reader.Data("../output/run" + str(idNum) + "/", nRead=4000)
	print(do.data[1].shape)

#	for i in range(do.nr):
#		print(i, do.r[i], do.tOrbit[i], do.data[20][0,i], do.data[19][0,i])

	#do.makeAnimFrames()

	# initial conditions
	do.icPlots()
	
	# 4 pannel space-time plots
	do.multiStPlot1(10, vLineCoords1=[50, 50+float(sys.argv[1])/2.0])

	# all time scales on the same plot (at t=0)
	n=0
	plt.loglog(do.r, do.tOrbit, label="tOrbit")
	plt.loglog(do.r, do.data[19][n], label="tAdv")
	plt.loglog(do.r, do.data[20][n], label="tDiff")
	plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
	plt.ylim(3.e-4,1.e6)
	plt.ylabel("timescale (years)")
	plt.xlabel("r (AU)")
	plt.savefig(do.savePath + "timeScales.png", bbox_inches='tight')
	plt.clf()

	# sigma multi profile
	do.multiProfile(15, 10, logOption=1, save="png", spacing="log")
	# alpha multi profile
	do.multiProfile(0, 10, logOption=1, save="png", spacing="log")
	# h/r multi profile
	do.multiProfile(1, 10, logOption=1, ymin=1.e-2, ymax=1.e0, save="png", spacing="log")
	# Tc^4 multi profile
	do.multiProfile(2, 10, logOption=1, ymin=1.e10, ymax=1.e19, save="png", spacing="lin")
	# nu multi profile
	do.multiProfile(7, 10, logOption=1, save="png", spacing="log")
	# mdot multi profile
	do.multiProfile(9, 10, logOption=1, save="png", spacing="log")
	# bz multi profile
	do.multiProfile(12, 10, logOption=0, ymin=-1.e1, ymax=1.e1, save="png", spacing="log")
	do.multiProfile(12, 10, logOption=1, ymin=-1.e1, ymax=1.e1, save="png", spacing="lin")
	do.multiProfile(12, 10, logOption=1, ymin=1.e-3, ymax=1.e1, save="png")
	# beta multi profile
	do.multiProfile(16, 10, logOption=1, save="png")

	# brs multi profile
	#do.multiProfile(6, 10, logOption=0, save="pdf", ymin=-1.0, ymax=1.0)
	#do.multiProfile(6, 10, logOption=1, save="pdf", ymin=1.e-3, ymax=1.e3)
	# psi multi profile
	do.multiProfile(14, 10, logOption=1, spacing="log", save="pdf")
	# brs/bz multi profile
	do.multiProfile(17, 10, logOption=0, save="pdf", ymin=-2.0, ymax=2.0)




	#do.profile(0, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	#do.profile(1, 0.01, vLineCoords=[10.0, 100.0], logOption=1, ymin=1.e-2, ymax=1.e0, save='png')
	#do.profile(9, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	#do.profile(12, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	#do.profile(15, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')

	
	#tStart1 = 10; tEnd1 = 32*10;
	'''
	tStart1 = None; tEnd1 = None;
	dt = float(sys.argv[1])
	vLineCoords1 = np.zeros(int(1+(do.tmax-50.0)/dt))
	vLineCoords2 = np.zeros(int(1+(do.tmax-50.0)/dt))
	i=0	
	print(int(1+(do.tmax-50.0)/dt))
	while i < len(vLineCoords1):
		vLineCoords1[i] = 50.0 + i*dt
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
	
	
	print(do.lum[-1]/do.lum[-100])
	'''


'''
	lamMaxMax = 0.0; lamMaxMin = 1.e40; specMaxMax = 0.0; specMaxMin = 1.e40
	for r in [0.1, 1.0, 10.0]:
		time = 100.0		
		rIndex = do.getrindex(r)
		temp = np.power(do.data[3][do.gettindex(100.0),rIndex],0.25)
		spectrum = do.calcSingleSpectrum(temp)
		specMax = np.amax(spectrum)
		lamMax = do.lam[np.argmax(spectrum)]
		if lamMax > lamMaxMax: lamMaxMax = lamMax
		if lamMax < lamMaxMin: lamMaxMin = lamMax
		if specMax > specMaxMax: specMaxMax = specMax
		if specMax < specMaxMin: specMaxMin = specMax
		plt.loglog(do.lam, spectrum, label="r="+str(np.round(r,2))+" AU")
	plt.xlim(lamMaxMin*7.e-2, lamMaxMax*2.e2); 
	plt.ylim(specMaxMin*1.e-5, specMaxMax*2.e0);
	plt.xlabel("wavelength (code units)"); plt.ylabel("Specific Intensity (code units)"); plt.legend()
	plt.savefig(do.savePath + "spectra.png", bbox_inches='tight'); plt.clf();


	lamMaxMax = 0.0; lamMaxMin = 1.e40; specMaxMax = 0.0; specMaxMin = 1.e40
	sumSpectrum = np.zeros_like(do.lam)
	for r in do.r:
		time = 100.0		
		rIndex = do.getrindex(r)
		temp = np.power(do.data[3][do.gettindex(100.0),rIndex],0.25)
		spectrum = do.calcSingleSpectrum2(time, r)
		sumSpectrum += spectrum
		specMax = np.amax(spectrum)
		lamMax = do.lam[np.argmax(spectrum)]
		if lamMax > lamMaxMax: lamMaxMax = lamMax
		if lamMax < lamMaxMin: lamMaxMin = lamMax
		if specMax > specMaxMax: specMaxMax = specMax
		if specMax < specMaxMin: specMaxMin = specMax
		if rIndex%20==0 and rIndex>10 :
			plt.loglog(do.lam, spectrum, label="r="+str(np.round(r,2))+" AU")
	plt.loglog(do.lam, sumSpectrum, label="sum")
	plt.xlim(lamMaxMin*7.e-2, lamMaxMax*2.e2); 
	plt.ylim(specMaxMin*1.e-5, specMaxMax*1.e2);
	plt.xlabel("wavelength (code units)"); plt.ylabel("Integrated Intensity (code units)");
	plt.legend(bbox_to_anchor=(1.02, 1.0), loc=2, borderaxespad=0.)
	plt.savefig(do.savePath + "spectra2.png", bbox_inches='tight'); plt.clf();

	for tIndex in range(0, int(do.nt/2), int(do.nt/2/10)):	
		sumSpectrum = np.zeros_like(do.lam)
		time = do.t[tIndex]	
		for r in do.r:		
			rIndex = do.getrindex(r)
			spectrum = do.calcSingleSpectrum2(time, r)
			sumSpectrum += spectrum
		plt.loglog(do.lam, sumSpectrum, label="t="+str(np.round(time,2))+" years")
	plt.xlim(1.e-5, 1.e-1); 
	plt.ylim(1.e10, 1.e18);
	plt.xlabel("wavelength (code units)"); plt.ylabel("Integrated Intensity (code units)");
	plt.legend(bbox_to_anchor=(1.02, 1.0), loc=2, borderaxespad=0.)
	plt.savefig(do.savePath + "spectra3.png", bbox_inches='tight'); plt.clf();

	for tIndex in range(0, int(do.nt), int(do.nt/100)):	
		sumSpectrum = np.zeros_like(do.lam)
		time = do.t[tIndex]	
		for r in do.r:		
			rIndex = do.getrindex(r)
			spectrum = do.calcSingleSpectrum2(time, r)
			sumSpectrum += spectrum
		plt.loglog(do.lam, sumSpectrum)
		plt.title("t="+str(np.round(time,2))+" years")
		plt.xlim(1.e-5, 1.e-1); 
		plt.ylim(1.e10, 1.e15);
		plt.xlabel("wavelength (code units)"); plt.ylabel("Integrated Intensity (code units)");
		plt.savefig(do.savePath + "spectra4_" + str(tIndex) + ".png", bbox_inches='tight'); plt.clf();
'''






'''
if str(sys.argv[2])=="save": 
	t0 = time.time()
	do.makeMultiAnim(save="yes", timeCutFactor=int(sys.argv[3]), lengthInSeconds=int(60))
	print(time.time() - t0)
if str(sys.argv[2])=="show": do.makeMultiAnim(show="yes")
'''


























