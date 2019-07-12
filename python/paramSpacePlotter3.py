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
def tDiff_rdz(alpha, mDot, mStar):
	alphaFac = np.power(alpha/1.e-2,-11./9.)
	mStarFac = np.power(mStar/1.0,1./3.)
	mDotFac  = np.power(mDot/1.e-7,4./9.)
	return 201. * alphaFac * mStarFac * mDotFac
tCycle = 20
################################################################################
colors = [(1,0,0,1),(0.8,0,0.2,1),(0.5,0,0.5,1),(0.2,0,0.8,1),(0,0,1,1)]
tCycle = 20.0
mDotArr = np.asarray([np.power(10,float(i)) for i in np.arange(-10,-6.9,0.1)])
mStarList = [np.power(10.,i) for i in np.arange(-1, 1.1, 0.5)]
labelList = mStarList
################################################################################
alpha = 0.03
n_mStar = 0
for mStar in mStarList:
	color = colors[n_mStar]
	tDiffArr = tDiff_rdz(alpha, mDotArr, mStar)
	plt.loglog(mDotArr, tDiffArr, color=color,
			   label=r"$M_{star}=$"+str(np.round(mStarList[n_mStar],1)))
	n_mStar+=1
plt.legend()
plt.title(r"$\alpha=$"+str(alpha))
plt.ylim(1,200)
plt.xlim(1.e-10,1.e-7)
plt.axhline(y=tCycle, color='k', linestyle='--')
plt.ylabel(r"$t_{cycle}$" + " (years)")
plt.xlabel(r"$\dot{m}$" + " " + r"$(M_\odot$" + "/year)")
plt.tight_layout()
plt.savefig("./paramSpace3.png", bbox_inches='tight')
plt.close('all')
################################################################################
alpha = 1.e-2
n_mStar = 0
for mStar in mStarList:
	color = colors[n_mStar]
	tDiffArr = tDiff_rdz(alpha, mDotArr, mStar)
	plt.loglog(mDotArr, tDiffArr, color=color,
			   label=r"$M_{star}=$"+str(np.round(mStarList[n_mStar],1)))
	n_mStar+=1
plt.legend()
plt.title(r"$\alpha=$"+str(alpha))
plt.ylim(1,200)
plt.xlim(1.e-10,1.e-7)
plt.axhline(y=tCycle, color='k', linestyle='--')
plt.ylabel(r"$t_{cycle}$" + " (years)")
plt.xlabel(r"$\dot{m}$" + " " + r"$(M_\odot$" + "/year)")
plt.tight_layout()
plt.savefig("./paramSpace4.png", bbox_inches='tight')
plt.close('all')
################################################################################















#
