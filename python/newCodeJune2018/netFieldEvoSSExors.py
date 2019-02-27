#!/usr/bin/python
import numpy as np
import os
import sys
import resource
import time
from shutil import copyfile
import functionLib as lib

# process input file
#inPath     = "../../input/"
#inFileName = "exors"+str(sys.argv[1])

inpArgs    = sys.argv[1:]
inp        = lib.InputFile(int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]))





tmax    = inp.tWait + inp.nCycles * inp.tCycle + inp.tRelax
tmax1   = inp.tWait + inp.nCycles * inp.tCycle
dtOut   = tmax/float(inp.nOut)

Tdz4    = 1.e15 #new
prandtl = 1.0

reportCutFactor = 500;
writeCutFactor  = 5000;

innerAdvBc  = 0
innerDiffBc = 1
outerAdvBc  = 1 
outerDiffBc = 1





tempGrid = np.loadtxt("../../fmatrix/outGrid_" + str(inp.gridId) + ".csv", delimiter=',')
nrOgGrid = len(tempGrid)
#for i in range(len(tempGrid)):
#	print(str(i) + ", " + str(tempGrid[i]))

riMinReal = (np.abs(tempGrid-inp.rMin)).argmin()
rMinReal  = tempGrid[riMinReal]
riInReal  = (np.abs(tempGrid-inp.rIn)).argmin()
rInReal   = tempGrid[riInReal]
riOutReal = (np.abs(tempGrid-inp.rOut)).argmin()
rOutReal  = tempGrid[riOutReal]
riMaxReal = (np.abs(tempGrid-inp.rMax)).argmin()
rMaxReal  = tempGrid[riMaxReal]

tempGrid2 = tempGrid[riMinReal:riMaxReal+1]

nrReal    = len(tempGrid2)


riDz1    = (np.abs(tempGrid2-inp.rDz1)).argmin()
#if riDz1<20: riDz1+=inp.nSmooth+1
riDz2    = (np.abs(tempGrid2-inp.rDz2)).argmin()

#print(riMinReal, rMinReal)
#print(riInReal,  rInReal)
#print(riOutReal, rOutReal)
#print(riMaxReal, rMaxReal)
#print(nrReal)
#print(riDz1)
#print(riDz2)

for i in range(len(tempGrid2)):
	if tempGrid2[i]<rInReal:	 riBuffer1=i+1
	if tempGrid2[i]<rOutReal:  riBuffer2=i+1			

print(tempGrid2.shape)



mu          = 1.0
mp          = 1.0
kr0         = 1.0
kb          = 1.0
rootG       = (2.0*3.14159)
rootGM      = rootG * np.sqrt(inp.mStar)
littleSigma = 1.0
rStarAu     = inp.rStar * 0.00465

savePath = "../../output/run"+str(inp.runId)+"/"
if not os.path.exists(savePath): os.makedirs(savePath)
inpArgFile = open(savePath+"params.txt", "w")
for thing in inpArgs:
	inpArgFile.write(thing)
	inpArgFile.write('\n')
inpArgFile.close()









# IC function for bz 
def getPsi0(sg):
	bz0 = inp.bInitScale*np.power(sg.r, inp.bz0index)
	return getPsiFromBz(sg,0.0,bz0) 

def getBz0(sg):
	bz0 = inp.bInitScale*np.power(sg.r, inp.bz0index)
	return bz0

def getStellarDipole(sg):
	b                 = inp.bStar*np.power(sg.r/rStarAu  , -3.0)
	bFlux             = b * sg.r * sg.dr
	bFluxTot          = np.sum(bFlux)
	bFluxPerBuffer    = bFluxTot / float(riBuffer1)
	bAdd              = np.zeros_like(b)
	bAdd[0:riBuffer1] = bFluxPerBuffer / (sg.r[0:riBuffer1]*sg.dr[0:riBuffer1])
	return bAdd

def getStellarDipoleRampUp(sg, teff):
	#bStar1 = inp.bStar * (teff/(tmax-inp.tWait))
	bStar1 = np.sign(inp.bStar)*min(np.absolute(inp.bStar), np.absolute(inp.bStar) * (teff/(0.5*inp.tmax-inp.tWait)) )
	b = bStar1*np.power(sg.r/rStarAu , -3.0)
	return b

# alpha calculation 
alphaSmoothingMatrix = np.zeros([nrReal,nrReal])
if inp.nSmooth > 0:
	for i in range(nrReal):
		if i<riBuffer1:
			alphaSmoothingMatrix[i,riBuffer1:riBuffer1+inp.nSmooth]=1.0
		elif riBuffer1<=i and i<riBuffer2:
			alphaSmoothingMatrix[i,max(i-inp.nSmooth,riBuffer1):min(i+inp.nSmooth+1,riBuffer2)]=1.0
		elif riBuffer2<=i and i<nrReal:
			alphaSmoothingMatrix[i,riBuffer2-inp.nSmooth:riBuffer2] = 1.0
else:
	for i in range(nrReal):
		alphaSmoothingMatrix[i,i]=1.0
for i in range(nrReal):
	alphaSmoothingMatrix[i]=alphaSmoothingMatrix[i]/np.sum(alphaSmoothingMatrix[i])
def getDzEdgeIndex(dg):
	return (np.abs(dg.Tc4-Tdz4)).argmin()
def getAlpha(sg, s):
	# alpha if active for all r
	alphaActive = np.clip(11.0*np.power(np.abs(s.beta), -0.53),inp.alphaMinAz,inp.alphaMaxAz)	
	# alpha if dead for all r	
	alphaDead = inp.alphaDz*np.ones_like(sg.r)
	# 1 if active, 0 if dead for all r
	activeOnes = np.zeros_like(sg.r)
	# 1 if dead, 0 if active for all r	
	deadOnes = np.zeros_like(sg.r)
	# inside riDz1 is always active
	activeOnes[0:riDz1]   = 1.0
	# Hall DZ can be active or dead
	activeOnes[riDz1:riDz2] = s.beta[riDz1:riDz2]>0.0
	deadOnes[riDz1:riDz2]   = s.beta[riDz1:riDz2]<0.0
	# if inbetween, dead if previously dead, active if previously active
	# active 
	activeOnes[riDz1:riDz2] = np.logical_or( 
															np.logical_and(
																np.logical_and(
																	s.bz[riDz1:riDz2] > -inp.threshFactor*np.abs(sg.bz0[riDz1:riDz2]),
																	s.bz[riDz1:riDz2] <  inp.threshFactor*np.abs(sg.bz0[riDz1:riDz2])															
															),
																s.alphaRawPrev[riDz1:riDz2] == inp.alphaMinAz
															),
															s.bz[riDz1:riDz2] > inp.threshFactor*np.abs(sg.bz0[riDz1:riDz2])
														)
		
	# dead
	deadOnes[riDz1:riDz2] = np.logical_or( 
															np.logical_and(
																np.logical_and(
																	s.bz[riDz1:riDz2] > -inp.threshFactor*np.abs(sg.bz0[riDz1:riDz2]),
																	s.bz[riDz1:riDz2] <  inp.threshFactor*np.abs(sg.bz0[riDz1:riDz2])															
																),
																s.alphaRawPrev[riDz1:riDz2] == inp.alphaDz
															),
															s.bz[riDz1:riDz2] < -inp.threshFactor*np.abs(sg.bz0[riDz1:riDz2])
														)
		
	#print(deadOnes[riDz1:riDz2])
	# outside riDz2 is always active
	activeOnes[riDz2:nrReal]  = 1.0
	#print(s.alphaRawPrev)
	#print(activeOnes)
	#print(deadOnes)
	# assign correct alpha for all r
	alphaRaw = activeOnes*alphaActive + deadOnes*alphaDead
	# smooth with smoothing matrix 
	alphaSmooth = np.dot(alphaSmoothingMatrix, alphaRaw)	
	return alphaRaw, alphaSmooth 

counter=0
# driving sigma and/or bz
def getAddToSig(sg, dg, s):
	addToSig = np.zeros_like(sg.r)
	return addToSig
def getAddToBz(sg, dg, s):
	global counter
	addToBz = np.zeros_like(sg.r)
	#addToBz[0:riBuffer1] = driveAmp
	teff = sg.t[-1]-inp.tWait
	tModCycle = teff%inp.tCycle
	if sg.t[-1]>inp.tWait and sg.t[-1]<tmax1:
		if inp.rampUpOption != 1:	
			if   (0.00*inp.tCycle) < tModCycle < (0.25*inp.tCycle): addToBz =  getStellarDipole(sg)*(dg.dt/(inp.tCycle/4.0))
			elif (0.25*inp.tCycle) < tModCycle < (0.50*inp.tCycle): addToBz =  np.zeros_like(getStellarDipole(sg))
			elif (0.50*inp.tCycle) < tModCycle < (0.75*inp.tCycle):	addToBz = -getStellarDipole(sg)*(dg.dt/(inp.tCycle/4.0))
			else                                          : addToBz =  np.zeros_like(getStellarDipole(sg))
		if inp.rampUpOption == 1:
			if   (0.00*inp.tCycle) < tModCycle < (0.25*inp.tCycle): addToBz =  getStellarDipoleRampUp(sg, teff)*(dg.dt/(inp.tCycle/4.0))
			elif (0.25*inp.tCycle) < tModCycle < (0.50*inp.tCycle): addToBz =  np.zeros_like(getStellarDipoleRampUp(sg, teff))
			elif (0.50*inp.tCycle) < tModCycle < (0.75*inp.tCycle):	addToBz = -getStellarDipoleRampUp(sg, teff)*(dg.dt/(inp.tCycle/4.0))
			else                                          : addToBz =  np.zeros_like(getStellarDipoleRampUp(sg, teff))
	#if teff < 10.0*inp.tCycle/2.0 :	return addToBz/2.0
	if teff < inp.tCycle/2.0 and inp.firstCycleFactor>0 :	return addToBz*inp.firstCycleFactor
	#else: return addToBz
	#if counter%100==0: print((getStellarDipoleRampUp(sg, teff))[10]);
	counter+=1;
	return addToBz
	

# functions to advance psi and sigma
def getkPsi(sg, dg, s):
	return -(sg.r*dg.vAdv*s.bz + sg.r*dg.vDiff*s.brs)  
def getkSig(sg, dg, s):
	returnArray = np.zeros_like(sg.r)
	nuSi   = dg.nu * sg.x * s.sig
	nuSip1 = np.roll(nuSi, -1)
	nuSim1 = np.roll(nuSi,  1)
	returnArray[0]    = np.square(1.0/(sg.x[0]    * sg.dx[0]))    * (nuSip1[0]                   - 2.0*nuSi[0])  
	returnArray[1:-1] = np.square(1.0/(sg.x[1:-1] * sg.dx[1:-1])) * (nuSip1[1:-1] + nuSim1[1:-1] - 2.0*nuSi[1:-1])
	returnArray[-1]   = np.square(1.0/(sg.x[-1]   * sg.dx[-1]))   * (inp.sigBcFactor*((sg.x[-1]+sg.dx[-1])/sg.x[-1])*nuSim1[-1] + nuSim1[-1] - 2.0*nuSi[-1])
	return returnArray
def rkGetNextState(sg, dg, s1):
	s2 = State(sg, dg, s1.psi + s1.kPsi*(dg.dt/2.0), s1.sig + s1.kSig*(dg.dt/2.0))
	s3 = State(sg, dg, s1.psi + s2.kPsi*(dg.dt/2.0), s1.sig + s2.kSig*(dg.dt/2.0))
	s4 = State(sg, dg, s1.psi + s3.kPsi*(dg.dt/1.0), s1.sig + s3.kSig*(dg.dt/1.0))
	kPsi = (s1.kPsi + 2.0*s2.kPsi + 2.0*s3.kPsi + s4.kPsi)/6.0;
	kSig = (s1.kSig + 2.0*s2.kSig + 2.0*s3.kSig + s4.kSig)/6.0;
	addToPsi = getPsiFromBz(sg, dg, getAddToBz(sg, dg, s1))
	addToSig = getAddToSig(sg, dg, s1)
	return State(sg, dg, s1.psi + kPsi*dg.dt + addToPsi, s1.sig + kSig*dg.dt + addToSig)

# helper stuff
def spaceDeriv(sg, a):
	aPlus=np.roll(a,1)
	aMinus=np.roll(a,-1)
	da = (aMinus-aPlus)/2.0
	da[0] = a[1] - a[0]
	da[sg.nr-1] = (a[sg.nr-1] - a[sg.nr-2]) 
	return da/sg.dr
def xIn(r):
	return (r-rMinReal)/(rInReal-rMinReal)
def xOut(r):
	return (r-rMaxReal)/(rOutReal-rMaxReal)
def smooth(x):
	return (3.0*x*x)/(1.0+2.0*x*x*x)
summingHelperMatrix = np.zeros([nrReal,nrReal])
for i in range(nrReal):
	summingHelperMatrix[i,0:i]=1.0
def getPsiFromBz(sg, dg, bz):
	temp = sg.r * sg.dr * bz
	psi = np.dot(summingHelperMatrix, temp)
	return psi
def getrindex(sg, r1):
	return (np.abs(sg.r-r1)).argmin()
def getTimeStep(sg, dg, s):
	dtAdv  = np.amin(np.abs(sg.dr/dg.vAdv))
	dtDiff = np.amin(np.abs(2.0*sg.dr/(3.14159*dg.vDiff)))
	dt = min(dtAdv, dtDiff, 0.1)
	return dt*inp.courantNo

	

# static grid object
class StaticGrid:
	def __init__(self, idNum):
		self.rMin  = rMinReal
		self.rIn   = rInReal
		self.rOut  = rOutReal
		self.rMax  = rMaxReal
		self.rTemp = np.loadtxt("../../fmatrix/outGrid_" + str(idNum) + ".csv", delimiter=',')
		self.r     = self.rTemp[riMinReal:riMaxReal+1]
		self.rootr = np.power(self.r, 0.5)
		self.x     = self.rootr
		self.dx    = np.zeros_like(self.r)
		self.nr    = len(self.r)		
		self.dr    = np.zeros_like(self.r)
		self.Omega = rootGM*np.power(self.r,-1.5)
		self.Omega2= np.power(self.Omega,2.0)
		self.fInvMatrixTemp = np.loadtxt("../../fmatrix/fmatrixInvFlat_" + str(idNum) +".csv", delimiter=',')
		self.fInvMatrixTemp = np.reshape(self.fInvMatrixTemp, [nrOgGrid, nrOgGrid])
		self.fInvMatrix     = self.fInvMatrixTemp[riMinReal:riMaxReal+1, riMinReal:riMaxReal+1]
		self.t = [0.0];
		self.bz0 = getBz0(self)
		for i in range(self.nr):
			if i==0: self.dr[i]=self.r[i+1]-self.r[i]
			elif 0<i and i<self.nr-1: self.dr[i]=(self.r[i+1]-self.r[i-1])/2.0
			elif i==self.nr-1: self.dr[i]=(self.r[i]-self.r[i-1])
		for i in range(self.nr):
			if i==0: self.dx[i]=self.x[i+1]-self.x[i]
			elif 0<i and i<self.nr-1: self.dx[i]=(self.x[i+1]-self.x[i-1])/2.0
			elif i==self.nr-1: self.dx[i]=(self.x[i]-self.x[i-1])
		self.applyAdvBc = np.zeros_like(self.r)		
		for i in range(self.nr):
			if self.r[i]<=self.rIn:
				if innerAdvBc==-1: self.applyAdvBc[i] = 0.0
				if innerAdvBc== 0: self.applyAdvBc[i] = smooth(xIn(self.r[i]))
				if innerAdvBc== 1: self.applyAdvBc[i] = 1.0 
			elif self.rIn<self.r[i] and self.r[i]<self.rOut: 
				self.applyAdvBc[i] = 1.0
			elif self.rOut<=self.r[i]:
				if outerAdvBc==-1: self.applyAdvBc[i] = 0.0
				if outerAdvBc== 0: self.applyAdvBc[i] = smooth(xOut(self.r[i]))
				if outerAdvBc== 1: self.applyAdvBc[i] = 1.0 
		self.applyDiffBc = np.zeros_like(self.r)		
		for i in range(self.nr):
			if self.r[i]<=self.rIn:
				if innerDiffBc==-1: self.applyDiffBc[i] = 0.0
				if innerDiffBc== 0: self.applyDiffBc[i] = smooth(xIn(self.r[i]))
				if innerDiffBc== 1: self.applyDiffBc[i] = 1.0 
			elif self.rIn<self.r[i] and self.r[i]<self.rOut: 
				self.applyDiffBc[i] = 1.0
			elif self.rOut<=self.r[i]:
				if outerDiffBc==-1: self.applyDiffBc[i] = 0.0
				if outerDiffBc== 0: self.applyDiffBc[i] = smooth(xOut(self.r[i]))
				if outerDiffBc== 1: self.applyDiffBc[i] = 1.0 

# state object
class State:
	def __init__(self, sg, dg, psi0, sig0, initialize=0):
		if initialize==1:
			self.psi = getPsi0(sg)
			dPsidr = spaceDeriv(sg, self.psi)
			self.bz = dPsidr/sg.r
			self.brs = np.dot(sg.fInvMatrix, self.psi)
			self.sig = np.zeros_like(sg.r)			
			self.kPsi = np.zeros_like(sg.r)
			self.kSig = np.zeros_like(sg.r)
			self.beta = 1.e10*np.ones_like(sg.r)
			self.alphaRawPrev = np.zeros_like(sg.r)
		elif initialize==0:
			self.psi = psi0
			self.sig = sig0
			dPsidr = spaceDeriv(sg, self.psi)
			self.bz = dPsidr/sg.r
			self.brs = np.dot(sg.fInvMatrix, self.psi)
			self.kPsi = getkPsi(sg, dg, self)
			self.kSig = getkSig(sg, dg, self)
			self.beta = np.sign(self.bz)*(dg.rho*np.power(dg.cs,2.0))/np.square(self.bz)
			self.alphaRawPrev = dg.alphaRaw

# dynamic grid object
class DynamicGrid:
	def __init__(self, sg, s, initialize=0):
		if initialize==1:
			if inp.bInitScale > 0.0: 
				self.alphaRaw = np.ones_like(sg.r)*inp.alphaMinAz
				self.alphaSmooth = np.ones_like(sg.r)*inp.alphaMinAz
			else:
				self.alphaRaw = np.ones_like(sg.r)*inp.alphaMinAz
				self.alphaRaw[riDz1:riDz2] = inp.alphaDz
				self.alphaSmooth = np.dot(alphaSmoothingMatrix, self.alphaRaw)
			self.mdot=np.ones_like(sg.r)*inp.mdot0
			#print(self.alphaSmooth)
			#print(self.alphaRaw)
			# self consistently solved quantities (9)
			self.Tc4    = np.power(np.power(1.5,0.2)*np.power(inp.mdot0,0.4)*np.power(mp,0.2)*np.power(kr0,0.2)*np.power(mu,0.2)*np.power(sg.Omega,0.6)*np.power(2.0,-1)*np.power(kb,-0.2)*np.power(3.14159,-0.4)*np.power(self.alphaSmooth,-0.2)*np.power(littleSigma,-0.2),4.0)
			self.Tdisk4 = np.power(np.power(2.0,0.75)*np.power(kb,0.25)*np.power(3.14159,0.25)*np.power(self.Tc4,5.0/16.0)*np.power(self.alphaSmooth,0.25)*np.power(inp.mdot0,-0.25)*np.power(mp,-0.25)*np.power(kr0,-0.25)*np.power(mu,-0.25)*np.power(sg.Omega,-0.25),4.0)
			self.cs     = np.sqrt(inp.mdot0*self.Tdisk4*kr0*sg.Omega/(self.Tc4*self.alphaSmooth)/(8.0*3.14159))
			self.rho    = (inp.mdot0*sg.Omega2)/(np.sqrt(18.0)*np.power(3.14159,1.5)*self.alphaSmooth*np.power(self.cs,3.0))
			self.h      = 2*self.cs/sg.Omega ############# THIS FACTOR OF 2
			self.kR     = np.sqrt(32.0/3.14159)*self.Tc4/(3.0*self.h*self.Tdisk4*self.rho)
			self.nu     = self.cs*self.h*self.alphaSmooth
			self.tau    = 4.0*self.Tc4/(3.0*self.Tdisk4)
			s.sig       = self.h*np.sqrt(2.0*3.14159)*self.rho
			# other intermediate helpful quantities	
			self.vAdv  = -1.5*(self.nu/sg.r)*sg.applyAdvBc
			self.vDiff = ((self.nu/sg.r)/(prandtl*(self.h/sg.r)))*sg.applyDiffBc
			s.kPsi      = getkPsi(sg, self, s)
			s.kSig      = getkSig(sg, self, s)
			s.alphaRawPrev = self.alphaRaw
			#plt.loglog(sg.r, self.alphaSmooth); plt.xlabel("r"); plt.ylabel("alpha"); plt.show(); plt.clf();
			#plt.loglog(sg.r, self.nu); plt.xlabel("r"); plt.ylabel("nu"); plt.show(); plt.clf();
			#plt.loglog(sg.r, self.mdot); plt.xlabel("r"); plt.ylabel("mdot"); plt.show(); plt.clf();
			#plt.loglog(sg.r, s.sig); plt.xlabel("r"); plt.ylabel("sigma"); plt.show(); plt.clf();
		elif initialize==0:
			self.alphaRaw, self.alphaSmooth = getAlpha(sg, s)
			# self consistently solved quantities (9)
			self.Tc4    = np.power((3.0/4.0)*np.power(kb,1.0/3.0)*np.power(self.alphaSmooth,1.0/3.0)*np.power(kr0,1.0/3.0)*np.power(s.sig,2.0/3.0)*np.power(sg.Omega,1.0/3.0)*np.power(mp,-1.0/3.0)*np.power(mu,-1.0/3.0)*np.power(littleSigma,-1.0/3.0) ,4.0)	
			self.Tdisk4 = np.power(np.power(2.0,0.75)*np.power(self.Tc4,0.25)*np.power(3.0,-0.25)*np.power(kr0,-0.25)*np.power(s.sig,-0.25),4.0)
			self.cs     = np.sqrt( (kb*np.power(self.Tc4,0.25)) / (mp*mu) )
			self.rho    = 9.0*self.cs*self.alphaSmooth*np.power(s.sig,2.0)*sg.Omega2 / (np.sqrt(128*3.14159)*self.Tdisk4*littleSigma)
			self.h      = 2*s.sig/(np.sqrt(2.0*3.14159)*self.rho) ############# THIS FACTOR OF 2
			self.kR     = 8.0*self.Tc4 / ( 3.0*self.Tdisk4*s.sig )
			self.nu     = self.cs*self.h*self.alphaSmooth
			self.tau    = 4.0*self.Tc4/(3.0*self.Tdisk4)
			self.mdot   = 3.0*3.14159*self.nu*s.sig
			# other intermediate helpful quantities	
			#self.vAdv  = -1.5*(self.nu/sg.r)*sg.applyAdvBc
			takeDerivOfThis = self.nu * s.sig * np.sqrt(sg.r)
			self.vAdv  = -(3.0/(s.sig*np.sqrt(sg.r))) * spaceDeriv(sg, takeDerivOfThis) * sg.applyAdvBc
			self.vDiff = ((self.nu/sg.r)/(prandtl*(self.h/sg.r)))*sg.applyDiffBc
			#plt.loglog(sg.r, self.alphaSmooth); plt.xlabel("r"); plt.ylabel("alpha"); plt.show(); plt.clf();
			#plt.loglog(sg.r, self.nu); plt.xlabel("r"); plt.ylabel("nu"); plt.show(); plt.clf();
			#plt.loglog(sg.r, self.mdot); plt.xlabel("r"); plt.ylabel("mdot"); plt.show(); plt.clf();
			#plt.loglog(sg.r, s.sig); plt.xlabel("r"); plt.ylabel("sigma"); plt.show(); plt.clf();
		self.dt = getTimeStep(sg, self, s)
		
# timer object
class Timer:
	def __init__(self, labels):
		self.nTimers = len(labels)
		self.count =      [0   for i in range(self.nTimers)]
		self.timeSpent =  [0.0 for i in range(self.nTimers)]
		self.startTime = [[]  for i in range(self.nTimers)]
		self.labels = labels
	def startTimer(self, n):
		self.count[n] = self.count[n]+1
		self.startTime[n].append(time.time())
	def endTimer(self, n):
		self.timeSpent[n] = self.timeSpent[n] + (time.time()-self.startTime[n][self.count[n]-1])
	def checkTimer(self, n):
		return time.time()-self.startTime[n][self.count[n]-1]
	def printInfo(self):
		for n in range(self.nTimers):
			if self.count[n]>1:
				sys.stdout.write(self.labels[n] + ": " + str(1000*self.timeSpent[n]/self.count[n]) + " ms per call  " + str(self.count[n]) + " calls" + "\n")
			elif self.count[n]==1:
				sys.stdout.write(self.labels[n] + ": " + str(1000*self.timeSpent[n]/self.count[n]) + " ms per call  " + str(self.count[n]) + " call" + "\n" )

# report function
def report(timer, n, t, tmax):
	msPerCycle = 1000*timer.checkTimer(0)/(max(n,1))
	msPerCycleRecent = 1000*(time.time()-timer.startTime[1][max(n-reportCutFactor,0)])/(reportCutFactor)
	timeRemaining = (timer.checkTimer(0)/(max(t,0.001))) * (tmax-t)	
	sys.stdout.write( "" "\n")
	sys.stdout.write( "runID = " + str(inp.runId) + "\n")
	sys.stdout.write( "cycle = " + str(n) + "   t = " + str(round(t,1)) + "/" + str(tmax) + "   " + str(round((t/tmax)*100,2)) + " % done" + "\n")
	sys.stdout.write( str(round(msPerCycle,3)) + " ms per cycle total" + "\n")
	sys.stdout.write( str(round(msPerCycleRecent,3)) + " ms per cycle recently " + "\n")    
	sys.stdout.write( str(int(timeRemaining/3600)) + ":" + str(int(timeRemaining%3600/60)) + ":" + str(int((timeRemaining%60))) + " remaining" + "\n") 
	sys.stdout.write("h/r at 1 AU = " + str(dg.h[(np.abs(sg.r-1.0)).argmin()]) )
	sys.stdout.flush()

# function to write to file
def writeToFile(sg, dgOut, sOut, nOutCurrent, timer):
	# manipulate into saveable output and save
	sgSaveArray = np.zeros([3, sg.nr])
	col=0;     sgSaveArray[col]=sg.r
	col=col+1; sgSaveArray[col]=sg.dr
	col=col+1; sgSaveArray[col]=sg.Omega
	np.save(savePath+"sgrid.npy", sgSaveArray)
	timer.startTimer(4)
	dgSaveArray = np.zeros([12, nOutCurrent, sg.nr])
	for n in range(len(dgOut)):
		col=0;     dgSaveArray[col,n]=dgOut[n].alphaSmooth
		col=col+1; dgSaveArray[col,n]=dgOut[n].h/sg.r
		col=col+1; dgSaveArray[col,n]=dgOut[n].Tc4
		col=col+1; dgSaveArray[col,n]=dgOut[n].Tdisk4
		col=col+1; dgSaveArray[col,n]=dgOut[n].cs
		col=col+1; dgSaveArray[col,n]=dgOut[n].rho
		col=col+1; dgSaveArray[col,n]=dgOut[n].kR
		col=col+1; dgSaveArray[col,n]=dgOut[n].nu
		col=col+1; dgSaveArray[col,n]=dgOut[n].tau
		col=col+1; dgSaveArray[col,n]=dgOut[n].mdot
		col=col+1; dgSaveArray[col,n]=dgOut[n].vAdv
		col=col+1; dgSaveArray[col,n]=dgOut[n].vDiff
	np.save(savePath+"dgrid.npy", dgSaveArray)
	timer.endTimer(4)
	timer.startTimer(5)
	sSaveArray = np.zeros([5, nOutCurrent, sg.nr])
	for n in range(len(dgOut)):
		col=0;     sSaveArray[col,n]=sOut[n].bz
		col=col+1; sSaveArray[col,n]=sOut[n].brs
		col=col+1; sSaveArray[col,n]=sOut[n].psi
		col=col+1; sSaveArray[col,n]=sOut[n].sig
		col=col+1; sSaveArray[col,n]=sOut[n].beta
	np.save(savePath+"state.npy", sSaveArray)
	timer.endTimer(5)
	np.save(savePath+"time.npy", np.asarray(tOut))




sg = StaticGrid(inp.gridId)

dgOut = []; sOut = []; tOut=[];

s = State(sg, 0.0, 0.0, 0.0, initialize=1)
dg = DynamicGrid(sg, s, initialize=1)

# integration loop
timer = Timer(["total", "cycles", "dgrid", "state", "saving dgrid", "saving state"])
n=0; tCheck = 0.0; nOutCurrent=0; 
timer.startTimer(0)
while tCheck<tmax:
	timer.startTimer(1)
	if (n)%reportCutFactor==0 and tCheck>0:
		report(timer, n, tCheck, tmax)
	if tCheck>=nOutCurrent*dtOut:
		sOut.append(s); dgOut.append(dg); tOut.append(sg.t[n]); nOutCurrent=nOutCurrent+1;
	if ( (n)%writeCutFactor==0 and n>0 ) :
		writeToFile(sg, dgOut, sOut, nOutCurrent, timer)
	n=n+1
	timer.startTimer(2)	
	dg = DynamicGrid(sg, s)
	timer.endTimer(2)
	timer.startTimer(3)
	s = rkGetNextState(sg, dg, s)
	timer.endTimer(3)
	sg.t.append(sg.t[n-1]+dg.dt)
	tCheck = sg.t[n]
	timer.endTimer(1)
timer.endTimer(0)
nOut=len(tOut)

writeToFile(sg, dgOut, sOut, nOutCurrent, timer)


timer.printInfo()








































