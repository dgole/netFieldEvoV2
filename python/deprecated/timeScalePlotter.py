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


r = np.logspace(np.log10(0.05), np.log10(6.0), num=400)
rootGM = 2*3.14
tOmega = r / (rootGM / np.sqrt(r))
tOrbit = tOmega * (2*3.14) 

hor = 0.1
alphaA = np.zeros_like(r)
alphaD = np.zeros_like(r)
for i in range(len(r)):
	if   i < 160: alphaD[i] = 4.e-2
	else        : alphaD[i] = 1.e-3
	alphaA[i] = 4.e-2

tAdvA  = tOmega * np.power(hor**2 * alphaA, -1)
tAdvD  = tOmega * np.power(hor**2 * alphaD, -1)
tDiffA = tOmega * np.power(hor    * alphaA, -1)
tDiffD = tOmega * np.power(hor    * alphaD, -1)

plt.loglog(r, tOrbit, label="orbit", color='k')
plt.loglog(r, tDiffA, label="diffusion, active", color='r')
plt.loglog(r, tAdvA,  label="advection, active", color='b')

plt.loglog(r, tDiffD, label="diffusion, dead", color='r', linestyle='--')
plt.loglog(r, tAdvD,  label="advection, dead", color='b', linestyle='--')


plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
plt.xlabel("r (AU)")
plt.ylabel("timescale (years)")
plt.savefig('./timeScales.png', bbox_inches='tight', dpi = 1000)
plt.clf()





















