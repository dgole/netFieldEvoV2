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
font = {'family' : 'sans-serif', 'weight' : 'normal', 'size'   : 16}
fontLabels = {'family' : 'sans-serif', 'weight' : 'normal', 'size'   : 20}
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
#header = [r"$\alpha$", r"$h/r$", r"$T_c^4$", r"$T_disk^4$", r"$c_s$", r"$\rho$", r"$\kappa_R$", r"$\nu$", r"$\tau$", r"$mdot$", r"$v_{adv}$", r"$v_{diff}$", r"$\beta$", r"$B_z$", r"$B_{rs}$", r"$\psi$", r"$\Sigma$"

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
tiMin = do.gettindex(t1)
tiMax = do.gettindex(t2)

nRows=4
nCols=2
figSize=(10,15)
widthRatios=[1,0.02]
tmin=t1
tmax=t2
fig   = plt.figure(figsize=figSize)
gs    = gridspec.GridSpec(nrows=nRows, ncols=2, height_ratios=np.ones(nRows), width_ratios=widthRatios)
ax    = []
cax   = []


for n in range(nRows):
	if n==0: ax.append(fig.add_subplot(gs[n, 0]))
	else   : ax.append(fig.add_subplot(gs[n, 0], sharex=ax[0]))
for n in range(nRows-1):
	cax.append(fig.add_subplot(gs[n, 1]))

extent = [tmin, tmax, np.log10(do.rmin), np.log10(do.rmax)]
aspect = 0.3*(tmax-tmin)/(np.log10(do.rmax)-np.log10(do.rmin))

axNum = 0
col   = 12
im     = ax[axNum].imshow(
						np.transpose(np.fliplr(do.data[col][tiMin:tiMax])),
						extent=extent,
						aspect=aspect,
						cmap=plt.get_cmap('coolwarm'),
						norm=colors.SymLogNorm(linthresh=0.1, linscale=1.0, vmin=-1.0, vmax=1.0)
						)
ax[axNum].set_ylabel('log(r) (AU)')
ax[axNum].set_title(do.header[col] + " (G)")
ax[axNum].axhline(np.log10(do.inp.rDz1), color='k', linestyle='--')
#ax[axNum].set_xlabel("Time (yr)")
fig.colorbar(im, cax=cax[axNum], orientation='vertical')

axNum = 1
col   = 9
im     = ax[axNum].imshow(
						np.transpose(np.fliplr(do.data[col][tiMin:tiMax])),
						extent=extent,
						aspect=aspect,
						cmap=plt.get_cmap('viridis'),
						norm=colors.LogNorm()
						)
ax[axNum].set_ylabel('log(r) (AU)')
ax[axNum].set_title(do.header[col])
ax[axNum].axhline(np.log10(do.inp.rDz1), color='k', linestyle='--')
#ax[axNum].set_xlabel("Time (yr)")
fig.colorbar(im, cax=cax[axNum], orientation='vertical')

axNum = 2
col   = 15
im     = ax[axNum].imshow(
						np.transpose(np.fliplr(do.data[col][tiMin:tiMax])),
						extent=extent,
						aspect=aspect,
						cmap=plt.get_cmap('viridis'),
						norm=colors.LogNorm()
						)
ax[axNum].set_ylabel('log(r) (AU)')
ax[axNum].set_title(do.header[col])
ax[axNum].axhline(np.log10(do.inp.rDz1), color='k', linestyle='--')
#ax[axNum].set_xlabel("Time (yr)")
fig.colorbar(im, cax=cax[axNum], orientation='vertical')

axNum = 3
smoothLum = np.zeros_like(do.lum)
for i in range(1, smoothLum.shape[0]-1):
	smoothLum[i] = (do.lum[i+1] + do.lum[i] + do.lum[i-1])/3.0

ax[axNum].semilogy(do.t[tiMin:tiMax], smoothLum[tiMin:tiMax])
ax[axNum].set_xlim(tmin, tmax)
ax[axNum].set_xlabel("Time (yr)")
ax[axNum].set_ylabel('Normalized Disk Luminosity')




plt.savefig(do.savePath + "fullST.png", bbox_inches='tight')



























#
