#!/usr/bin/python
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import os
import math
import sys
from matplotlib.backends.backend_pdf import PdfPages
import time
################################################################################
def rDipole(r0, lam):
	return r0 * np.square(np.cos(lam))
def toCart(lam, r):
	x = r * np.cos(lam)
	y = r * np.sin(lam)
	return x, y
def adjustX(x, xMax):
	xMax2 = xMax * 0.96
	xAdj = (xMax2*2.0)/(1+np.exp(-(x/(xMax2)))) - xMax2
	return xAdj
################################################################################
xMax = 0.2
xIn  = np.arange(0,1,0.05)
xAdj = adjustX(xIn, xMax)
plt.plot(xIn, xAdj)
plt.xlim(0,1)
plt.ylim(0,xMax)
plt.savefig('../diagram2.png'); plt.clf()
################################################################################
#fig = plt.figure(figsize=(4,13), dpi=80)
#ax = []
#ax.append(plt.subplot2grid((3, 1), (0, 0), rowspan=1))
#ax.append(plt.subplot2grid((3, 1), (1, 0), rowspan=1))
#ax.append(plt.subplot2grid((3, 1), (2, 0), rowspan=1))
################################################################################
fig = plt.figure(figsize=(13,4), dpi=100)
ax = []
ax.append(plt.subplot2grid((1, 3), (0, 0), rowspan=1))
ax.append(plt.subplot2grid((1, 3), (0, 1), rowspan=1))
ax.append(plt.subplot2grid((1, 3), (0, 2), rowspan=1))
################################################################################
# plot original stellar dipole
lam    = np.arange(0, 2*3.2, 0.01)
r0arr  = np.power(np.arange(4.0,9,0.5), 3.0)
r0arr /= (np.amax(r0arr)*1.0)
xMax   = 0.2
for r0 in r0arr:
	r    = rDipole(r0, lam)
	x, y = toCart(lam, r)
	ax[0].plot(x, y, color=(0,0,1,1), linewidth=1)
################################################################################
# plot adjusted stellar dipole
for r0 in r0arr:
	r    = rDipole(r0, lam)
	x, y = toCart(lam, r)
	x    = adjustX(x, xMax)
	ax[2].plot(x, y, color=(0,0,1,1), linewidth=1)
################################################################################
# plot open field lines
pltLim = xMax*2.2
yArr   = np.arange(-pltLim*1.5, pltLim*1.5, 0.1)
xArr   = np.arange(xMax*0.8, pltLim, 0.05)
for x in xArr:
	ax[1].plot(x*np.ones_like(yArr), yArr, color='r', linewidth=1)
	ax[1].plot(-x*np.ones_like(yArr), yArr, color='r', linewidth=1)
xArr   = [0.20, 0.22, 0.25, 0.29, 0.34, 0.4]
for x in xArr:
	ax[2].plot(x*np.ones_like(yArr), yArr, color='r', linewidth=1)
	ax[2].plot(-x*np.ones_like(yArr), yArr, color='r', linewidth=1)
################################################################################
# plot central star
ax[0].plot(0.0, 0.0, 'ko', markersize=20 )
ax[1].plot(0.0, 0.0, 'ko', markersize=20 )
ax[2].plot(0.0, 0.0, 'ko', markersize=20 )
################################################################################
# plot disk
x = np.arange(xMax*0.75, pltLim, 0.01)
y = 0.15 * x
ax[1].fill_between( x, y1=-y, y2=y, color=(0,0,0,0.2))
ax[1].fill_between(-x, y1=-y, y2=y, color=(0,0,0,0.2))
ax[2].fill_between( x, y1=-y, y2=y, color=(0,0,0,0.2))
ax[2].fill_between(-x, y1=-y, y2=y, color=(0,0,0,0.2))
################################################################################
for thisAx in ax:
	thisAx.set_aspect('equal')
	thisAx.set_xlim(-pltLim, pltLim)
	thisAx.set_ylim(-pltLim, pltLim)
	thisAx.get_xaxis().set_visible(False)
	thisAx.get_yaxis().set_visible(False)
	thisAx.set_facecolor((0,0,0,0.05))
################################################################################
ax[0].set_title('Stellar Dipole: Closed Lines')
ax[1].set_title('Disk Net-Field: Open Lines')
ax[2].set_title('Star + Disk Fields')
#ax[2].annotate('closed stellar lines \n become open disk lines',
 			   #xy=(xMax*0.9, 0.01),
			   #xytext=(-0.17, -pltLim*0.99))
#ax[2].annotate("",
            #xy=(xMax*0.93, 0.02),
            #xytext=(0.03, -pltLim*0.83),
            #arrowprops=dict(facecolor='black', shrink=0.05, width=2))
#ax[2].annotate("",
            #xy=(-xMax*0.93, 0.02),
            #xytext=(-0.03, -pltLim*0.83),
            #arrowprops=dict(facecolor='black', shrink=0.05, width=2))
################################################################################
fig.tight_layout()
fig.savefig('./diagram_horizontal.png'); plt.clf()
################################################################################
