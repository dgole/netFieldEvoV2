#!/usr/bin/python
import numpy as np
import os
import sys
import resource
import time
import matplotlib.pyplot as plt

au = 1.5e13
G  = 6.67e-8
sigma=5.67e-5
Msun=2e33
Rsun=0.00465*au
year=3.15e7
days=24.0*60.0*60.0

class InputFile:
	def __init__(self, runId, mStar, bStar, mdot0, tCycle, nCycles):
		self.runId = runId
		self.mStar = mStar
		self.bStar=bStar
		self.tCycle=tCycle
		self.nCycles=nCycles
		self.mdot0  = mdot0
		self.gridId=8001
		self.tWait=200.0
		self.tRelax=0
		self.nOut=8000
		self.alphaMaxAz = 1.e-1
		self.alphaMinAz = 1.e-1
		self.alphaDz    = 1.e-1 ################################
		#self.alphaDz    = 3.e-3
		self.bInitScale = -1.e-2  ###############################
		#self.bInitScale = -1.e-3  ###############################
		#self.bInitScale = -1.e-6  ###############################
		self.threshFactor = 0.0
		self.nSmooth      = 5
		self.rampUp       = 0
		self.firstCycleFactor = 0.5
		self.courantNo    = 0.4
		self.bz0index     = -0.5
		self.sig0option   = 0
		self.sigOindex    = 0.0
		self.rDz2         = 100000
		self.sigBcFactor  = 0.9935
		self.rampUpOption = 0
		self.pStar        = 7.0
		if self.mStar==1.0:
			self.rStar = 2.0
		elif self.mStar==0.9:
			self.rStar = 1.75
		elif self.mStar==0.8:
			self.rStar = 1.5
		elif self.mStar==0.7:
			self.rStar = 1.5
		elif self.mStar==0.6:
			self.rStar = 1.8
		elif self.mStar==0.5:
			self.rStar = 1.8
		elif self.mStar==0.4:
			self.rStar = 1.6
		elif self.mStar==0.3:
			self.rStar = 2.0
		elif self.mStar==0.2:
			self.rStar = 1.5
		elif self.mStar==0.1:
			self.rStar = 1.0


		self.rMs  = 7.93e-13 * np.power(np.square(self.bStar) * np.power(self.rStar*Rsun,6) * np.power(self.mdot0*Msun/year,-1) * np.power(self.mStar*Msun,-0.5), 2.0/7.0)
		self.rCo  = np.power(G*(self.mStar*Msun)*np.square(self.pStar*days)/(4*3.1415*3.1415), 1.0/3.0) / au
		self.rIn = min(self.rMs, self.rCo)

		r     = np.linspace(self.rIn,10,num=1.e4)
		td    = np.power((3*G*self.mStar*Msun*self.mdot0*Msun/year)/(8.0*3.1415*sigma*np.power(r*au,3))*(1-np.sqrt(self.rStar*Rsun/(r*au))),0.25)
		kr    = 3.0
		sig   = (self.mStar)*1700*np.power(r,-3.0/2.0)
		tau   = 0.5 * sig * kr
		tc    = np.power(np.power(td,4.0)*0.75*tau,0.25) / 1.0
		index = (np.abs(tc-1000)).argmin()
		self.rDz1 = r[index]

		self.rMin = self.rIn *(1.0/1.25)
		self.rOut = min(self.rIn * 100.0, 5.0/1.25)
		self.rMax = self.rOut*(1.25)

		#self.rDz1 = self.rMax #####################################################



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
