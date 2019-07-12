#!/usr/bin/python
import numpy as np
import os
import sys
import resource
import time

# process input file111
inp = np.asarray(np.genfromtxt("../input/largeDisk"+str(sys.argv[1])+".txt", dtype=str))
print(inp)
for i in range(inp.shape[0]):
	word = inp[i][0]
	number = inp[i][1]
	if   word == 'runId'        : runId        = int(number)
	elif word == 'gridId'       : gridId       = int(number)
	elif word == 'nr'           : nr           = int(number)
	elif word == 'tmax'         : tmax         = float(number)
	elif word == 'rMin'         : rMin         = float(number)
	elif word == 'rMax'         : rMax         = float(number)
	elif word == 'rIn'          : rIn          = float(number)
	elif word == 'rOut'         : rOut         = float(number)
	elif word == 'rYearNorm'    : rYearNorm    = float(number)
	elif word == 'nOut'         : nOut         = float(number)
	elif word == 'courantNo'    : courantNo    = float(number)
	elif word == 'mdot0'        : mdot0        = float(number)
	elif word == 'bInitScale'   : bInitScale   = float(number)
	elif word == 'bz0index'     : bz0index     = float(number)
	elif word == 'alpha0'       : alpha0       = float(number)
	elif word == 'sigBcFactor'  : sigBcFactor  = float(number)
	elif word == 'prandtl'      : prandtl      = float(number)
	#elif word == '' :  = number

Tdz4 = 1.e15 #new
reportCutFactor = 500;
writeCutFactor = 5000;
innerAdvBc  = 0
innerDiffBc = 1
outerAdvBc  = 0 
outerDiffBc = 0
dtOut = tmax/float(nOut)
tempGrid = np.loadtxt("../fmatrix/outGrid_" + str(gridId) + ".csv", delimiter=',')
for i in range(nr):
	print(str(i) + ", " + str(tempGrid[i]))
for i in range(nr):
	if tempGrid[i]:
		if tempGrid[i]<rIn:
			riBuffer1=i+1
		if tempGrid[i]<rOut:
			riBuffer2=i+1			
mu          = 1.0
mp          = 1.0
kr0         = 1.0
kb          = 1.0
rootGM      = (2.0*3.14159)*np.power(rYearNorm,1.5)
yearTime    = (2.0*3.14159)*np.power(rYearNorm,1.5)/rootGM
print(yearTime)
littleSigma = 1.0

savePath = "../outputLargeDisk/run"+str(runId)+"/"
if not os.path.exists(savePath): os.makedirs(savePath)

# IC function for bz 
def getPsi0(sg):
	bz0 = bInitScale*np.power(sg.r, bz0index)
	return getPsiFromBz(sg,0.0,bz0) 

def getBz0(sg):
	bz0 = bInitScale*np.power(sg.r, bz0index)
	return bz0

# alpha calculation 
def getAlpha(sg, s):
	alpha = np.ones_like(sg.r)*alpha0
	return alpha 
	
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
	returnArray[-1]   = np.square(1.0/(sg.x[-1]   * sg.dx[-1]))   * (sigBcFactor*((sg.x[-1]+sg.dx[-1])/sg.x[-1])*nuSim1[-1] + nuSim1[-1] - 2.0*nuSi[-1])
	return returnArray
def rkGetNextState(sg, dg, s1):
	s2 = State(sg, dg, s1.psi + s1.kPsi*(dg.dt/2.0), s1.sig + s1.kSig*(dg.dt/2.0))
	s3 = State(sg, dg, s1.psi + s2.kPsi*(dg.dt/2.0), s1.sig + s2.kSig*(dg.dt/2.0))
	s4 = State(sg, dg, s1.psi + s3.kPsi*(dg.dt/1.0), s1.sig + s3.kSig*(dg.dt/1.0))
	kPsi = (s1.kPsi + 2.0*s2.kPsi + 2.0*s3.kPsi + s4.kPsi)/6.0;
	kSig = (s1.kSig + 2.0*s2.kSig + 2.0*s3.kSig + s4.kSig)/6.0;
	return State(sg, dg, s1.psi + kPsi*dg.dt, s1.sig + kSig*dg.dt)

# helper stuff
def spaceDeriv(sg, a):
	aPlus=np.roll(a,1)
	aMinus=np.roll(a,-1)
	da = (aMinus-aPlus)/2.0
	da[0] = a[1] - a[0]
	da[sg.nr-1] = (a[sg.nr-1] - a[sg.nr-2]) 
	return da/sg.dr
def xIn(r):
	return (r-rMin)/(rIn-rMin)
def xOut(r):
	return (r-rMax)/(rOut-rMax)
def smooth(x):
	return (3.0*x*x)/(1.0+2.0*x*x*x)
summingHelperMatrix = np.zeros([nr,nr])
for i in range(nr):
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
	dt = min(dtAdv, dtDiff)
	return dt*courantNo

	

# static grid object
class StaticGrid:
	def __init__(self, idNum):
		self.r = np.loadtxt("../fmatrix/outGrid_" + str(idNum) + ".csv", delimiter=',')
		self.rootr = np.power(self.r, 0.5)
		self.x = self.rootr
		self.dx = np.zeros_like(self.r)
		self.nr = len(self.r)		
		self.dr = np.zeros_like(self.r)
		self.rMin = rMin
		self.rIn = rIn
		self.rOut = rOut
		self.rMax = rMax
		self.Omega = rootGM*np.power(self.r,-1.5)
		self.Omega2 = np.power(self.Omega,2.0)
		self.fInvMatrix = np.loadtxt("../fmatrix/fmatrixInvFlat_" + str(idNum) +".csv", delimiter=',')
		self.fInvMatrix = np.reshape(self.fInvMatrix, [self.nr, self.nr])
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
			if self.r[i]<self.rIn:
				if innerAdvBc==-1: self.applyAdvBc[i] = 0.0
				if innerAdvBc== 0: self.applyAdvBc[i] = smooth(xIn(self.r[i]))
				if innerAdvBc== 1: self.applyAdvBc[i] = 1.0 
			elif self.rIn<self.r[i] and self.r[i]<self.rOut: 
				self.applyAdvBc[i] = 1.0
			elif self.rOut<self.r[i]:
				if outerAdvBc==-1: self.applyAdvBc[i] = 0.0
				if outerAdvBc== 0: self.applyAdvBc[i] = smooth(xOut(self.r[i]))
				if outerAdvBc== 1: self.applyAdvBc[i] = 1.0 
		self.applyDiffBc = np.zeros_like(self.r)		
		for i in range(self.nr):
			if self.r[i]<self.rIn:
				if innerDiffBc==-1: self.applyDiffBc[i] = 0.0
				if innerDiffBc== 0: self.applyDiffBc[i] = smooth(xIn(self.r[i]))
				if innerDiffBc== 1: self.applyDiffBc[i] = 1.0 
			elif self.rIn<self.r[i] and self.r[i]<self.rOut: 
				self.applyDiffBc[i] = 1.0
			elif self.rOut<self.r[i]:
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
		elif initialize==0:
			self.psi = psi0
			self.sig = sig0
			dPsidr = spaceDeriv(sg, self.psi)
			self.bz = dPsidr/sg.r
			self.brs = np.dot(sg.fInvMatrix, self.psi)
			self.kPsi = getkPsi(sg, dg, self)
			self.kSig = getkSig(sg, dg, self)
			self.beta = np.sign(self.bz)*(dg.rho*np.power(dg.cs,2.0))/np.square(self.bz)

# dynamic grid object
class DynamicGrid:
	def __init__(self, sg, s, initialize=0):
		if initialize==1:
			self.alphaSmooth = getAlpha(sg,s)
			self.mdot=np.ones_like(sg.r)*mdot0
			#print(self.alphaSmooth)
			#print(self.alphaRaw)
			# self consistently solved quantities (9)
			self.Tc4    = np.power(np.power(1.5,0.2)*np.power(mdot0,0.4)*np.power(mp,0.2)*np.power(kr0,0.2)*np.power(mu,0.2)*np.power(sg.Omega,0.6)*np.power(2.0,-1)*np.power(kb,-0.2)*np.power(3.14159,-0.4)*np.power(self.alphaSmooth,-0.2)*np.power(littleSigma,-0.2),4.0)
			self.Tdisk4 = np.power(np.power(2.0,0.75)*np.power(kb,0.25)*np.power(3.14159,0.25)*np.power(self.Tc4,5.0/16.0)*np.power(self.alphaSmooth,0.25)*np.power(mdot0,-0.25)*np.power(mp,-0.25)*np.power(kr0,-0.25)*np.power(mu,-0.25)*np.power(sg.Omega,-0.25),4.0)
			self.cs     = np.sqrt(mdot0*self.Tdisk4*kr0*sg.Omega/(self.Tc4*self.alphaSmooth)/(8.0*3.14159))
			self.rho    = (mdot0*sg.Omega2)/(np.sqrt(18.0)*np.power(3.14159,1.5)*self.alphaSmooth*np.power(self.cs,3.0))
			self.h      = self.cs/sg.Omega
			self.kR     = np.sqrt(32.0/3.14159)*self.Tc4/(3.0*self.h*self.Tdisk4*self.rho)
			self.nu     = self.cs*self.h*self.alphaSmooth
			self.tau    = 4.0*self.Tc4/(3.0*self.Tdisk4)
			s.sig       = self.h*np.sqrt(2.0*3.14159)*self.rho
			# other intermediate helpful quantities	
			self.vAdv  = -1.5*(self.nu/sg.r)*sg.applyAdvBc
			self.vDiff = ((self.nu/sg.r)/(prandtl*(self.h/sg.r)))*sg.applyDiffBc
			s.kPsi      = getkPsi(sg, self, s)
			s.kSig      = getkSig(sg, self, s)
			#plt.loglog(sg.r, self.alphaSmooth); plt.xlabel("r"); plt.ylabel("alpha"); plt.show(); plt.clf();
			#plt.loglog(sg.r, self.nu); plt.xlabel("r"); plt.ylabel("nu"); plt.show(); plt.clf();
			#plt.loglog(sg.r, self.mdot); plt.xlabel("r"); plt.ylabel("mdot"); plt.show(); plt.clf();
			#plt.loglog(sg.r, s.sig); plt.xlabel("r"); plt.ylabel("sigma"); plt.show(); plt.clf();
		elif initialize==0:
			self.alphaSmooth = getAlpha(sg, s)
			# self consistently solved quantities (9)
			self.Tc4    = np.power((3.0/4.0)*np.power(kb,1.0/3.0)*np.power(self.alphaSmooth,1.0/3.0)*np.power(kr0,1.0/3.0)*np.power(s.sig,2.0/3.0)*np.power(sg.Omega,1.0/3.0)*np.power(mp,-1.0/3.0)*np.power(mu,-1.0/3.0)*np.power(littleSigma,-1.0/3.0) ,4.0)	
			self.Tdisk4 = np.power(np.power(2.0,0.75)*np.power(self.Tc4,0.25)*np.power(3.0,-0.25)*np.power(kr0,-0.25)*np.power(s.sig,-0.25),4.0)
			self.cs     = np.sqrt( (kb*np.power(self.Tc4,0.25)) / (mp*mu) )
			self.rho    = 9.0*self.cs*self.alphaSmooth*np.power(s.sig,2.0)*sg.Omega2 / (np.sqrt(128*3.14159)*self.Tdisk4*littleSigma)
			self.h      = s.sig/(np.sqrt(2.0*3.14159)*self.rho)
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
	sys.stdout.write( "runID = " + str(sys.argv[1]) + "\n")
	sys.stdout.write( "cycle = " + str(n) + "   t = " + str(round(t,1)) + "/" + str(tmax) + "   " + str(round((t/tmax)*100,2)) + " % done" + "\n")
	sys.stdout.write( str(round(msPerCycle,3)) + " ms per cycle total" + "\n")
	sys.stdout.write( str(round(msPerCycleRecent,3)) + " ms per cycle recently " + "\n")    
	sys.stdout.write( str(int(timeRemaining/3600)) + ":" + str(int(timeRemaining%3600/60)) + ":" + str(int((timeRemaining%60))) + " remaining" + "\n") 
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




sg = StaticGrid(gridId)

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
	#dg = DynamicGrid(sg, s)	
	dg = DynamicGrid(sg, s, initialize=1)
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








































