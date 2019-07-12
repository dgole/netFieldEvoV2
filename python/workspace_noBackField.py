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

#m.rcParams['text.usetex'] = True
#m.rcParams['text.latex.unicode'] = True
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

idNum = 2011

# make data object
do = reader.Data("../../output/run" + str(idNum) + "/")

# time scales
#reader.timeScalesPlot(do)
#plt.axhline(do.inp.tCycle, linestyle='--', color='k')
#plt.axhline(1.0, linestyle='--', color='k')
#plt.axvline(do.inp.rDz1, linestyle='--', color='k')
#plt.savefig(do.savePath + "timeScales.png", bbox_inches='tight'); plt.clf()

################################################################################

nrCut = 0

bzPreFactor = 1.e-3
col = 12
bzRms = bzPreFactor * np.sqrt(np.mean(np.square(do.data[col][do.nt-4000:,:]),axis=0))
plt.loglog(do.r[:do.nr-nrCut], bzRms[:do.nr-nrCut])
plt.xlabel('r (AU)')
plt.ylabel(do.header[col])

for preFactor in [1.5e-4]:
    model = preFactor*np.power(do.r[:do.nr-50],-3.5)
    plt.loglog(do.r[:do.nr-50], model, color='k', linestyle='--')

plt.savefig(do.savePath + "bz_rms_pl.png", bbox_inches='tight'); plt.clf()

################################################################################

beta1 = np.abs(do.data[16][5, do.getrindex(3.0)])
print(beta1)
normFactor = (1.e5/beta1) * 1.e6 # stellar field enhance
print(beta1*normFactor)
plotData = normFactor*np.power(np.mean(np.power(np.absolute(do.data[16][do.nt-4000:,:]),-1), axis=0),-1)
plt.loglog(do.r[:do.nr-nrCut], plotData[:do.nr-nrCut])
plt.xlabel('r (AU)')
plt.ylabel(do.header[16])

plt.savefig(do.savePath + "betaz_rms_pl.png", bbox_inches='tight'); plt.clf()


















#
