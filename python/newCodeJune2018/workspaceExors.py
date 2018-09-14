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


tWait = 5.0

idNum = int(sys.argv[1])
tFlip = float(sys.argv[1])

# make data object
do = reader.Data("../../output/run" + str(idNum) + "/")
print(do.data[1].shape)

col=1
reader.profile(do, col, 10)
plt.savefig(do.savePath + "profile_" + str(col) + "_.png", bbox_inches='tight')
plt.clf()

col=12
reader.st(do, col, logOption=2)
plt.savefig(do.savePath + "st_" + str(col) + "_.png", bbox_inches='tight')
plt.clf()



# initial conditions
#do.icPlots()

# all time scales on the same plot (at t=0)
'''
n1=0
n2=-1
print("tOrbit at 1 AU:")
print(do.tOrbit[do.getrindex(1.0)])
print(do.tOrbit[125])
plt.loglog(do.r, do.tOrbit, label="orbital", color='k')
plt.loglog(do.r, do.data[19][n1], label="advection, dead", color='b', linestyle='--')
plt.loglog(do.r, do.data[20][n1], label="diffusion, dead", color='b')
plt.loglog(do.r, do.data[20][n2]*10, label="advection, active", color='r', linestyle='--')
plt.loglog(do.r, do.data[20][n2], label="diffusion, active", color='r')
plt.xlim(7.e-2, 5.e0)
plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
#plt.legend(loc='best')
plt.ylim(1.e-2,1.5e5)
plt.ylabel("timescale (years)", fontdict=fontLabels)
plt.xlabel("r (AU)", fontdict=fontLabels)
plt.savefig(do.savePath + "timeScales.png", bbox_inches='tight')
plt.clf()
'''

#	for i in range(do.nr):
#		print(i, do.r[i], do.tOrbit[i], do.data[20][0,i], do.data[19][0,i])

	#do.makeAnimFrames()


	
# 4 pannel space-time plots
#coords = [i for i in range(50,int(do.tmax),int(tFlip))]
#coords1 = [tWait, tWait+tFlip/2]
#do.multiStPlot1(tFlip, vLinesOn=False, vLineCoords1=[tWait, tWait+5])
#do.multiStPlot1cbar(tFlip, vLinesOn=False, vLineCoords1=[tWait, tWait+5])

#ti1 = do.gettindex(do.tmax-tFlip*12.1); ti2 = do.gettindex(do.tmax-tFlip*0.0);
#do.multiStPlot1(tFlip, ti1=ti1, ti2=ti2, nameTag='_short')
#do.multiStPlot1cbar(tFlip, ti1=do.nt-1000, ti2=do.nt-1, nameTag='_short')
#do.multiStPlot3(tFlip, vLinesOn=False)
#do.multiStPlot3(tFlip, ti1=do.nt-1000, ti2=do.nt-1, nameTag='_short')
#nCycles = (do.tmax-tWait)/(2.0*tFlip)

#ti1 = do.gettindex(do.tmax-tFlip*7.1); ti2 = do.gettindex(do.tmax-tFlip*2.6);
#do.multiStPlot1(tFlip, ti1=ti1, ti2=ti2, nameTag='_shorter')
#do.multiStPlot1cbar(tFlip, ti1=ti1, ti2=ti2, nameTag='_shorter')
#do.multiStPlot3(tFlip, ti1=ti1, ti2=ti2, nameTag='_shorter')
#do.multiStPlot2(tFlip)


#do.timeEvo(12, 1.0, logOption=1, save='png')

# sigma multi profile
#do.multiProfile(15, 10, logOption=1, save="png", spacing="log")
# alpha multi profile
#do.multiProfile(0, 10, logOption=1, save="png", spacing="log")
# h/r multi profile
#do.multiProfile(1, 10, logOption=1, ymin=1.e-2, ymax=1.e0, save="png", spacing="log")
# Tc^4 multi profile
#do.multiProfile(2, 10, logOption=1, ymin=1.e10, ymax=1.e19, save="png", spacing="lin")
# nu multi profile
#do.multiProfile(7, 10, logOption=1, save="png", spacing="log")
# mdot multi profile
#do.multiProfile(9, 10, logOption=1, save="png", spacing="log")
# bz multi profile
#do.multiProfile(12, 10, logOption=0, ymin=-1.e1, ymax=1.e1, save="png", spacing="log")
#do.multiProfile(12, 10, logOption=1, ymin=-1.e1, ymax=1.e1, save="png", spacing="lin")
#do.multiProfile(12, 10, logOption=1, ymin=1.e-3, ymax=1.e1, save="png")
# beta multi profile
#do.multiProfile(16, 10, logOption=1, save="png")

# brs multi profile
#do.multiProfile(6, 10, logOption=0, save="pdf", ymin=-1.0, ymax=1.0)
#do.multiProfile(6, 10, logOption=1, save="pdf", ymin=1.e-3, ymax=1.e3)
# psi multi profile
#do.multiProfile(14, 10, logOption=1, spacing="log", save="png")
# brs/bz multi profile
#do.multiProfile(17, 10, logOption=0, save="png", ymin=-2.0, ymax=2.0)


#do.profile(10, 1000, logOption=0, save='png',ymin=-1.e-1, ymax=1.e-1)
#do.profile(10, 1000, logOption=0, save='png',ymin=-1.e-1, ymax=1.e-1)
#do.profile(10, 1005, logOption=0, save='png',ymin=-1.e-1, ymax=1.e-1)
#do.profile(10, 1010, logOption=0, save='png',ymin=-1.e-1, ymax=1.e-1)

#do.profile(10, 1000, logOption=1, save='png',ymin=1.e-5, ymax=1.e-1)
#do.profile(10, 1000, logOption=1, save='png',ymin=1.e-5, ymax=1.e-1)
#do.profile(10, 1005, logOption=1, save='png',ymin=1.e-5, ymax=1.e-1)
#do.profile(10, 1010, logOption=1, save='png',ymin=1.e-5, ymax=1.e-1)
#do.profile(1, 0.01, vLineCoords=[10.0, 100.0], logOption=1, ymin=1.e-2, ymax=1.e0, save='png')
#do.profile(9, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
#do.profile(12, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
#do.profile(15, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')

	
#tStart1 = 10; tEnd1 = 32*10;
	


# MAKE SINGLE PANE ST PLOT FOR BZ FOR FIGURE 3
#ti1 = 0	
#ti2 = do.nt-1
#thisNt = ti2-ti1
#aspect1 = 15.3 * (do.t[ti2]-do.t[ti1])/do.nr

#col = 12
#axNum = 0
#title = do.header[col]		
#plotData = do.data[col][ti1:ti2] 
#plotData=plotData/np.amax(np.abs(plotData))
#plt.imshow(np.transpose(np.fliplr(plotData)), 
#					extent=[do.t[ti1],do.t[ti2-1],np.log10(do.rmin),np.log10(do.rmax)], 
#					aspect=aspect1, 
#					cmap='coolwarm', 
#					norm=colors.SymLogNorm(linthresh=0.003, linscale=0.5, vmin=-1.0, vmax=1.0)
#					)
#plt.ylabel("log(r) (AU)")
#plt.axvline(x=50, linestyle='--', color='k')
#plt.axvline(x=55, linestyle='--', color='k')
#plt.xlabel("t (years)")
#plt.title(r"$B_z$")
#plt.savefig(do.path + "STjustBz" + ".png", bbox_inches='tight', dpi = 600)
#plt.colorbar()
#plt.savefig(do.path + "STjustBz_cbar" + ".png", bbox_inches='tight', dpi = 600)
#plt.clf()





















