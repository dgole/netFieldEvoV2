#!/usr/bin/python
from __future__ import unicode_literals
import matplotlib as m
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import math
from scipy import fftpack as fft
import sys
from matplotlib.backends.backend_pdf import PdfPages
import resource
import time
import numpy.random
import matplotlib.animation as animation
from matplotlib.pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import functionLib as lib

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

colorNo = 'red';
colorTiny = 'green'
colorFull = 'blue' 
marker='o';



def addMarker(mDot, tCycle, color, marker):
	mDot   = [mDot]; tCycle = [tCycle]
	plt.loglog(mDot, tCycle, marker=marker, linestyle='None', color=color)

plt.ylabel(r'$t_{cycle}$')
plt.xlabel(r'$\dot{M}$')

# 0.1 Msun, alphaActive=0.03, alphaDead = 0.001
addMarker(1.e-7,  6.e0, colorNo, marker)
addMarker(1.e-8,  6.e0, colorNo, marker)
addMarker(1.e-9,  6.e0, colorNo, marker)
addMarker(1.e-10, 6.e0, colorNo , marker)

addMarker(1.e-7,  2.e1, colorNo, marker)
addMarker(1.e-8,  2.e1, colorNo, marker)
addMarker(1.e-9,  2.e1, colorNo, marker)
addMarker(1.e-10, 2.e1, colorTiny, marker)

addMarker(1.e-7,  6.e1, colorNo, marker)
addMarker(1.e-8,  6.e1, colorTiny, marker)
addMarker(1.e-9,  6.e1, colorTiny, marker)
addMarker(1.e-10, 6.e1, colorFull, marker)

addMarker(1.e-7,  2.e2, colorTiny, marker)
addMarker(1.e-8,  2.e2, colorFull, marker)
addMarker(1.e-9,  2.e2, colorFull, marker)
addMarker(1.e-10, 2.e2, colorFull, marker)








plt.savefig('./combinedPlot.png', bbox_inches='tight')















