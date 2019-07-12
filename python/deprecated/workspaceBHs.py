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

#profile(self, col, time, logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None)
#addCol(self, funcName, label, *args, **kwargs)

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

# make data object
do = reader.Data("../output/run" + str(sys.argv[1]) + "/")

print do.nt
#print do.data[19][10,86-5]
#print do.data[19][10,168-5]
#print do.data[20][10,86-5]
#print do.data[20][10,168-5]

# all time scales on the same plot
n=0
plt.loglog(do.r, do.tOrbit, label="tOrbit")
plt.loglog(do.r, do.data[19][n], label="tAdv")
plt.loglog(do.r, do.data[20][n], label="tDiff")
plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
plt.ylim(1.e-4,1.e5)
plt.ylabel("timescale (code units)")
plt.savefig(do.pdfName, format='pdf', bbox_inches='tight'); plt.clf()
'''
ri=80;  plt.plot(do.t, do.data[12][...,ri], label=str(ri))
ri=82;  plt.plot(do.t, do.data[12][...,ri], label=str(ri))
ri=84; plt.plot(do.t, do.data[12][...,ri], label=str(ri))
ri=86; plt.plot(do.t, do.data[12][...,ri], label=str(ri))
ri=88; plt.plot(do.t, do.data[12][...,ri], label=str(ri))
ri=90; plt.plot(do.t, do.data[12][...,ri], label=str(ri))
ri=100; plt.plot(do.t, do.data[12][...,ri], label=str(ri))
ri=110; plt.plot(do.t, do.data[12][...,ri], label=str(ri))
plt.axhline(0.0)
plt.xlabel("time"); plt.ylabel(r"$B_z/B_{z,0}$"); plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
plt.savefig(do.pdfName, format='pdf', bbox_inches='tight'); plt.clf()
'''
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
do.multiProfile(12, 10, logOption=1, save="pdf")
# beta multi profile
do.multiProfile(16, 10, logOption=1, save="pdf")
# temp multi profile
do.multiProfile(2, 10, logOption=1, save="pdf")


# brs multi profile
#do.multiProfile(6, 10, logOption=0, save="pdf", ymin=-1.0, ymax=1.0)
#do.multiProfile(6, 10, logOption=1, save="pdf", ymin=1.e-3, ymax=1.e3)
# psi multi profile
do.multiProfile(14, 10, logOption=1, spacing="log", save="pdf")
# brs/bz multi profile
do.multiProfile(17, 10, logOption=0, save="pdf", ymin=-2.0, ymax=2.0)






do.pdfName.close()



if str(sys.argv[2])=="save": 
	t0 = time.time()
	do.makeMultiAnim(save="yes", timeCutFactor=2, lengthInSeconds=60)
	print time.time() - t0
if str(sys.argv[2])=="show": do.makeMultiAnim(show="yes")



























