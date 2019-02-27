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
import matplotlib.gridspec as gridspec
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


idNumList = [1000, 1001, 1002, 1003]
tDiff_rDz = []

for idNum in idNumList:
	# make data object
	do = reader.Data("../../output/run" + str(idNum) + "/")

	# time scales
	rDz  = do.inp.rDz1
	rDzi = do.getrindex(rDz) - do.inp.nSmooth
	ti   = do.gettindex(do.inp.tWait) 
	#plt.loglog(do.r, do.data[20][ti])
	#plt.axvline(do.r[rDzi])	
	#plt.show()
	tDiff_rDz.append()
	

	#plt.axvline(do.inp.rDz1, linestyle='--', color='k')
	#plt.savefig(do.savePath + "timeScales.png", bbox_inches='tight'); plt.clf()

	#print(do.data[0].shape)
	#print(do.data[1][10,100])
























