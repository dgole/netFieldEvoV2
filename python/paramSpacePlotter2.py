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
hor0   = 0.0223
mdot0  = 1.e-7
alpha = 1.e-1
################################################################################
def omega(r, mStar):
	rootG  = (2.0*3.14159)
	rootGM = rootG * np.sqrt(mStar)
	return rootGM * np.power(r, -1.5)

def hor(mdot):
	return hor0 * np.power(mdot/mdot0, 0.2)

def t_diff(r, mStar, mdot, alpha):
	hor1   = hor(mdot)
	omega1 = omega(r, mStar)
	return np.power(omega1*alpha*hor1,-1)
################################################################################
'''
fig = plt.figure(figsize=(12,8), dpi=100)
ax = []
ax.append(plt.subplot2grid((2, 2), (0, 0), rowspan=1))
ax.append(plt.subplot2grid((2, 2), (0, 1), rowspan=1))
ax.append(plt.subplot2grid((2, 2), (1, 0), rowspan=1))
ax.append(plt.subplot2grid((2, 2), (1, 1), rowspan=1))
'''
axNum = 0
colors = [(1,0,0,1), (0.8,0,0.3,1), (0.4,0,0.7,1), (0,0,1,1)]
tCycle = 40.0
rDzMatrix = [[0.35,0.21,0.12,0.07],
			 [0.57,0.34,0.20,0.12],
			 [0.72,0.42,0.25,0.15],
			 [0.98,0.58,0.35,0.204]]
mDotList = [1.e-7, 1.e-8, 1.e-9, 1.e-10]
mStarList = [0.1, 0.3, 0.5, 1.0]
r = np.arange(0.0,1.01, 0.001)
tDiffList = []
for mStar in mStarList:
	n_mDot = 0
	for mDot in mDotList:
		color = colors[n_mDot]
		rDz   = rDzMatrix[axNum][n_mDot]
		tDiff = [t_diff(r1, mStar, mDot, alpha) for r1 in r]
		tDiff = np.asarray(tDiff)
		#ax[axNum].loglog(r, tDiff, color=color, label="mDot="+str(mDot))
		#arg = np.argmin(np.absolute(tDiff-tCycle))
		#ax[axNum].axvline(x=r[arg], color=color, linestyle='--')
		#ax[axNum].axvline(x=rDz, color=color, linestyle='--')
		#ax[axNum].set_ylabel(r"$t$")
		#ax[axNum].set_xlabel(r"$r$")
		n_mDot+=1
		rIndex = np.argmin(np.absolute(r-rDz))
		tDiffList.append(tDiff[rIndex])
	#plt.loglog(mDotList,tDiffList); plt.show(); plt.clf();
	#ax[axNum].axhline(y=tCycle, color='k', linestyle='--')
	#ax[axNum].set_title(mStar)
	#ax[axNum].set_ylim(0,150)
	#ax[axNum].set_title("Mstar = " + str(mStar))
	#ax[axNum].legend()
	axNum+=1
#plt.tight_layout()
#plt.savefig("./paraSpace.png", bbox_inches='tight')
#plt.close('all')
################################################################################
axNum = 0
for n in range(4):
	print("####################################")
	print(mStarList[axNum], mDotList[n])
	print(tDiffList[n])
	print( np.log10(tDiffList[n]/tDiffList[n+1]) )







################################################################################
