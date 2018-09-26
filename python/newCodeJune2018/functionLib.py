#!/usr/bin/python
import numpy as np
import os
import sys
import resource
import time

class InputFile:
	def __init__(self, fileName):
		inp = np.asarray(np.genfromtxt(fileName, dtype=str))
		print(inp)
		for i in range(inp.shape[0]):
			word = inp[i][0]
			number = inp[i][1]
			if   word == 'runId'         : self.runId        = int(number)
			elif word == 'gridId'        : self.gridId       = int(number)
			elif word == 'nr'            : self.nr           = int(number)
			elif word == 'rMin'          : self.rMin         = float(number)
			elif word == 'rMax'          : self.rMax         = float(number)
			elif word == 'rIn'           : self.rIn          = float(number)
			elif word == 'rOut'          : self.rOut         = float(number)
			elif word == 'mStar'         : self.mStar        = float(number)
			elif word == 'tWait'         : self.tWait        = float(number)
			elif word == 'tCycle'        : self.tCycle       = float(number)
			elif word == 'tRelax'        : self.tRelax       = float(number)
			elif word == 'nCycles'       : self.nCycles      = float(number)
			elif word == 'nOut'          : self.nOut         = float(number)
			elif word == 'courantNo'     : self.courantNo    = float(number)
			elif word == 'mdot0'         : self.mdot0        = float(number)
			elif word == 'bInitScale'    : self.bInitScale   = float(number)
			elif word == 'bz0index'      : self.bz0index     = float(number)
			elif word == 'sig0option'    : self.sig0option   = int(number)
			elif word == 'sig0index'     : self.sig0index    = float(number)
			elif word == 'bStar'         : self.bStar        = float(number)
			elif word == 'rStar'         : self.rStar        = float(number)
			elif word == 'threshFactor'  : self.threshFactor = float(number)
			elif word == 'riDz1'         : self.riDz1        = int(number)
			elif word == 'riDz2'         : self.riDz2        = int(number) 
			elif word == 'alphaMaxAz'    : self.alphaMaxAz   = float(number)
			elif word == 'alphaMinAz'    : self.alphaMinAz   = float(number)
			elif word == 'alphaDz'       : self.alphaDz      = float(number)
			elif word == 'sigBcFactor'   : self.sigBcFactor  = float(number)
			elif word == 'nSmooth'       : self.nSmooth      = int(number)
			elif word == 'rampUp'        : self.rampUpOption = int(number)
			elif word == 'halfFirstCycle': self.halfFirstCycleOption = int(number)
			#elif word == '' :  = number





























