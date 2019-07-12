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

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}
m.rc('font', **font)

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
# 10 -vadv
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

def saveAndClear(figNum=0):
	plt.figure(figNum)
	plt.savefig(pdfName, format='pdf', bbox_inches='tight');
	plt.clf()


for idNum in sys.argv[1:]: 
	do = reader.Data("../../outputLargeDisk/run" + str(idNum) + "/", nRead=8000)
	savePath = "./plots/"+str(idNum)+"/"
	if not os.path.exists(savePath): os.makedirs(savePath)
	pdfName = PdfPages(savePath + "all.pdf")
	######################################################################################
	for col in range(0,22):
		x, y = reader.getProfileData(do, col, 0.1)
		index = np.log10(y[do.getrindex(10)] / y[do.getrindex(1)])
		print(do.header[col] + ' power law index: ' + str(index))
	######################################################################################
	reader.profilePlot(do, 19, 1, color='r', label='header')
	reader.profilePlot(do, 20, 1, color='b', label='header')
	plt.legend(loc='best');
	saveAndClear()
	reader.multiProfilePlot(do, 17, 10, logOption='semilogx')
	plt.ylim(-2,2)
	saveAndClear()
	reader.multiProfilePlot(do, 12, 10)
	saveAndClear()
	reader.multiProfilePlot(do, 14, 10)
	saveAndClear()
	reader.multiProfilePlot(do, 16, 10)
	saveAndClear()
	######################################################################################
	for col in range(0,21):
		reader.multiProfilePlot(do, col, 10)
		saveAndClear()

	pdfName.close()















