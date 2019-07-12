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
m.rcParams['text.latex.unicode'] = True
font = {'family' : 'sans-serif', 'weight' : 'normal', 'size'   : 20}
fontLabels = {'family' : 'sans-serif', 'weight' : 'normal', 'size'   : 24}
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

idNum = int(sys.argv[1])

# make data object
do = reader.Data("../../output/run" + str(idNum) + "/")

t1 = float(sys.argv[2])
t2 = float(sys.argv[3])
n1 = do.gettindex(t1)
n2 = do.gettindex(t2)
print(n1, n2)

smoothLum = np.zeros_like(do.lum)
for i in range(2, smoothLum.shape[0]-2):
	smoothLum[i] = np.mean(do.lum[i-2:i+3])


plt.figure(figsize=(10,5))
plt.semilogy(do.t[n1:n2], smoothLum[n1:n2], color='k')
plt.xlim(t1, t2)
plt.ylim(6.e-1, 6.e1)
plt.xlabel("time (yr)")
plt.ylabel("Norm. Disk Lum.")


if do.inp.mStar == 1.0: starStr = 'G'
elif do.inp.mStar == 0.8: starStr = 'K'
elif do.inp.mStar == 0.5: starStr = 'M1'
elif do.inp.mStar == 0.3: starStr = 'M4'
elif do.inp.mStar == 0.1: starStr = 'M6'

mdotStr = r"$\dot{m}=$" + str(do.inp.mdot0) + r"$M_\odot$" + "/yr"

tCycleStr = r"$t_{cycle}$=" + str(do.inp.tCycle)

plt.title(starStr + " star,   " + tCycleStr )

fileName = starStr + "_" + str(int(-np.log10(do.inp.mdot0))) + "_" + str(int(do.inp.tCycle))
print(fileName)
plt.savefig(do.savePath + "justLum.png", bbox_inches='tight')
plt.savefig("../../output/paperFigs/justLum_" + fileName +".png", bbox_inches='tight')





































#
