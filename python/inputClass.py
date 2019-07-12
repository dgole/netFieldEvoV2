#!/usr/bin/python
import numpy as np
import os
import sys
import resource
import time
import matplotlib.pyplot as plt

# define constants
au = 1.5e13
G  = 6.67e-8
sigma=5.67e-5
Msun=2e33
Rsun=0.00465*au
year=3.15e7
days=24.0*60.0*60.0

#runId  mStar  bStar   bInit   mdot0  FCF  tCycle nCycles
class InputFile:
	'''
	Holds all relevant model parameters in an object.
	Passed to both the model itself and the analysis module.
	'''
	def __init__(self, runId, mStar, bInit, mdot0, firstCycleFactor, tCycle, nCycles):
		self.runId = runId
		self.mStar = mStar
		self.tCycle=tCycle
		self.nCycles=nCycles
		self.mdot0  = mdot0
		self.bInitScale = bInit
		self.firstCycleFactor = firstCycleFactor
		self.gridId=8001
		self.tWait=50.0
		self.tRelax=0
		self.nOut=8000
		self.alphaMaxAz = 1.e-2
		self.alphaMinAz = 1.e-2
		self.alphaDz    = 3.e-4
		self.threshFactor = 0.0
		self.nSmooth      = 5
		self.rampUpOption = 0
		self.courantNo    = 0.1
		self.bz0index     = -0.5
		self.sig0option   = 0
		self.sigOindex    = 0.0
		self.rDz2         = 100000
		self.sigBcFactor  = 0.9935

		if self.mStar==1.0:
			self.rStar = 2.0
			self.pStar = 7.0
			self.bStar = -1.e3

		elif self.mStar==0.8:
			self.rStar = 1.6
			self.pStar = 12.0
			self.bStar = -1.e3

		elif self.mStar==0.5:
			self.rStar = 1.0
			self.pStar = 2.0
			self.bStar = -2.e3

		elif self.mStar==0.3:
			self.rStar = 0.6
			self.pStar = 2.0
			self.bStar = -2.e3

		elif self.mStar==0.1:
			self.rStar = 0.4
			self.pStar = 2.0
			self.bStar = -2.e3


		self.rMs  = 7.93e-13 * np.power(np.square(self.bStar) * np.power(self.rStar*Rsun,6) * np.power(self.mdot0*Msun/year,-1) * np.power(self.mStar*Msun,-0.5), 2.0/7.0)
		self.rCo  = np.power(G*(self.mStar*Msun)*np.square(self.pStar*days)/(4*3.1415*3.1415), 1.0/3.0) / au

		self.rIn = min(self.rMs, self.rCo)

		self.rDz1 = np.power(1000./767.,-10./9.)*np.power(self.alphaMaxAz/1.e-2,-2./9.)*np.power(self.mdot0/1.e-7,4./9.)*np.power(self.mStar/1.0,1./3.)
		if self.rDz1 < self.rIn: self.rDz1 = self.rIn

		self.rMin = self.rIn *(1.0/1.25)
		self.rOut = min(self.rIn * 100.0, 5.0/1.25)
		self.rMax = self.rOut*(1.25)

		tDiff0 = 308
		tDiffRdz = tDiff0 * np.power(1000./767.,-29./18.) * np.power(self.alphaMaxAz/1.e-2,-11./9.) * np.power(self.mStar/1.,1./3.) * np.power(self.mdot0/1.e-7,4./9.)

		print("mStar = " + str(self.mStar))
		print("rStar = " + str(self.rStar))
		print('')
		print('rMs  = ' + str(self.rMs))
		print('rCo  = ' + str(self.rCo))
		print('')
		print('rMin = ' + str(self.rMin))
		print('rIn  = ' + str(self.rIn))
		print('rDz  = ' + str(self.rDz1))
		print('rOut = ' + str(self.rOut))
		print('rMax = ' + str(self.rMax))
		print('')
		print('tDiff at rDz from PLs = ' + str(tDiffRdz))
		print('')





#
