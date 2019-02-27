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

idNum = int(sys.argv[1])

# make data object
do = reader.Data("../../output/run" + str(idNum) + "/")

# multiple ST plots
mp1 = reader.MultiPannel()
mp1.addImshowBz(0, do)
mp1.addImshowOther(1, do, 0)
mp1.addImshowOther(2, do, 15)
mp1.addLum(3, do)

plt.savefig(do.savePath + "multiPannel.png", bbox_inches='tight')


'''

im = ax[0].imshow(
						      np.transpose(np.fliplr(do.data[12])),
						      extent=extent,
				      		  aspect=aspect,
				    		  cmap=plt.get_cmap('coolwarm'),
						      norm=colors.SymLogNorm(linthresh=0.01, linscale=1.0, vmin=-10.0, vmax=10.0)
					        )
ax[0].set_ylabel('log(r) (AU)')
ax[0].set_title(do.header[12])
fig.colorbar(im, cax=cax[0], orientation='vertical')
#, ticks=[-1.e-1, -1.e-3, 0, 1.e-3, 1.e-1])

im = ax[1].imshow(
		       			  np.transpose(np.fliplr(do.data[0])),
		 				      extent=extent,
						      aspect=aspect,
						      cmap=plt.get_cmap('viridis'),
							  norm=colors.LogNorm()
						      )
ax[1].set_ylabel('log(r) (AU)')
ax[1].set_title(do.header[0])
fig.colorbar(im, cax=cax[1], orientation='vertical')

im = ax[2].imshow(
						      np.transpose(np.fliplr(do.data[15])),
					          extent=extent,
						      aspect=aspect,
						      cmap=plt.get_cmap('viridis'),
							  norm=colors.LogNorm()
						      )
ax[2].set_ylabel('log(r) (AU)')
ax[2].set_title(do.header[15])
fig.colorbar(im, cax=cax[2], orientation='vertical')


ax[3].semilogy(do.t, do.lum)
ax[3].set_xlim(do.t[0],do.t[-1])
ax[3].set_xlabel('t (years)')
ax[3].set_ylabel('Luminosity')


plt.savefig(do.savePath + "multiPannel.png", bbox_inches='tight')
'''
























