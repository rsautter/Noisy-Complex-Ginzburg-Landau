import numpy as np
from numpy.fft import fftn,ifftn,fftfreq
import math
import random
import itertools
import tqdm as tqdm
from collections import namedtuple
from scipy.signal import convolve2d
import cNoise
from scipy.interpolate import interp1d

class NCGL():
	'''
	NCGL - Noisy Complex Ginzburg-Landau
	
	Wrote by: Rubens Andreas Sautter (2022)
	
	Adapted from Aranson, et.al.(1997)
	https://arxiv.org/abs/patt-sol/9709005
	
	Adittional References:

	de Franciscis, d’Onofrio (2012) (Tsallis-Borland model)
	https://journals.aps.org/pre/abstract/10.1103/PhysRevE.86.021118
	
	
	
	Complex Ginzburg-Landau equation solved with Fourier pseudospectral methods, and integrated with RK45.
	
	A new method 

	'''

	def __init__(self, c1=1.0, c2=1.0,h=1.0, msize = 128, ic='r', sigma_r= 1.0, noiseType='multiplicative', noiseArgs=None,dim=2):
		'''
		Spatial parameters:
			ic = initial condition('r', 'g')
			h - grid spacing
			dim - dimension of the data (integer)

		GL parameters:
			c1 - diffusion parameter - (1+ib)Nabla A
			c2 - reaction parameter  - (1+ic)(|A|^2)A
			
		Noise Parameters:
			noiseSpeed - ]0,1[ - the speed (relative to the number of iterations) which the noise moves
			sigma_r - reactive noise 'strenght'
			noiseArgs - Colored noise parameters {'beta':2,std = 0.01}
		'''
		
		self.c1, self.c2 = c1,c2
		self.a0 = 0.75
		
		self.h = h
		self.ic = ic
		self.msize = msize
		self.dim = dim
		self.sigma_r = sigma_r
		self.noiseType = noiseType
		if noiseArgs is None:
			self.noiseArgs = {}
		else:
			self.noiseArgs = noiseArgs

	def __getRandom(self,n,dim):
		newShape = tuple([n for i in range(dim)])
		return np.random.rand(n**dim).reshape(newShape)
		
	def __getGaussian(self,n,dim):
		out = np.zeros(np.repeat(n,dim))
		c = n/2
		squareDists = np.sum((np.indices(out.shape)-c)**2,axis=0)
		return 2*(np.exp(-squareDists/(20*n))-0.5)
		
	def getInitialCondition(self):
		if self.ic=='r':
			self.a = self.a0*((self.__getRandom(self.msize,self.dim)-0.5)+1j*(self.__getRandom(self.msize,self.dim)-0.5))
		else:
			self.a = self.a0*(self.__getGaussian(self.msize,self.dim)+1j*self.__getGaussian(self.msize,self.dim))
			
		return np.array(self.a)
		
	
	def getChainedSingleReaction(self,a0=None,dt=0.1, nit=3000):
		'''
		Returns the iteration of a single amplitude (spatial part is ignored)
		
		The function integrates with rk4 method
		'''
		states = []
		delta = 1e-6*(np.random.rand()-0.5)
		if a0 is None:
			at = self.a0+delta
		else:
			at = a0
				
		for i in range(nit):
			states.append(at)
			
			t = i*dt
			k1 = self.reaction(at,t)
			k2 = self.reaction((at+dt*k1/2), t+dt/2.)
			k3 = self.reaction((at+dt*k2/2), t+dt/2.)
			k4 = self.reaction((at+dt*k3), t+dt)
			at = at + dt*(k1+2*k2+2*k3+k4)/6.
		return np.array(states)

	'''		
	def getNoisyLyapReaction(self,a0=None,beta=0,dt=0.1, nit=3000):

		lyap = []
		delta = 1e-6*(np.random.rand()-0.5)
		if a0 is None:
			at = self.a0+delta
		else:
			at = a0
		
		eta = cNoise.cNoise(beta=beta,shape=(nit+2,),std=1)
		eta2 = eta.copy()
		eta = np.gradient(eta)
		
		
		last = a0
		state = np.identity(3)
		lyapSpec=[]
		deltas = []
		for i in range(nit):
			a = np.real(at)
			b = np.imag(at)
			c = self.c2
			
			t = i*dt
			if self.noiseType == 'multiplicative':
				jac = [
					[1-3*(a**2)-b**2+2*a*b*c+self.sigma_r*eta[i], 	-2*a*b+(a**2)*c+3*(b**2)*c,		self.sigma_r*a-self.sigma_r*b],
					[-2*a*b-3*(a**2)*c-(b**2)*c, 	1-a**2-3*(b**2)-2*a*b*c+self.sigma_r*eta[i],		self.sigma_r*a+self.sigma_r*b],
					[0, 		0, 			np.exp(eta2[i]) ]
				]
				
				k1 = self.reaction(at,t)+self.sigma_r*at*self.__interpolate1D(eta,i)
				k11 = np.dot(jac,state)
				
				k2 = self.reaction((at+dt*k1/2), t+dt/2.)+self.sigma_r*(at+dt*k1/2)*self.__interpolate1D(eta,i+0.5)
				k22 = np.dot(jac,state+dt*k11/2)
				
				k3 = self.reaction((at+dt*k2/2), t+dt/2.)+self.sigma_r*(at+dt*k2/2)*self.__interpolate1D(eta,i+0.5)
				k33 = np.dot(jac,state+dt*k22/2)
				
				k4 = self.reaction((at+dt*k3), t+dt)+self.sigma_r*(at+dt*k3)*self.__interpolate1D(eta,i+1)
				k44 = np.dot(jac,state+dt*k33)
			else:
				jac = [
					[1-3*(a**2)-b**2+2*a*b*c, 	-2*a*b+(a**2)*c+3*(b**2)*c,		self.sigma_r],
					[-2*a*b-3*(a**2)*c-(b**2)*c, 	1-a**2-3*(b**2)-2*a*b*c,		0],
					[0, 		0, 				np.exp(eta2[i])  ]
				]
				k1 = self.reaction(at,t)+self.sigma_r*self.__interpolate1D(eta,i)				
				k11 = np.dot(jac,state)
				
				k2 = self.reaction((at+dt*k1/2), t+dt/2.)+self.sigma_r*self.__interpolate1D(eta,i+0.5)
				k22 = np.dot(jac,state+dt*k11/2)
				
				k3 = self.reaction((at+dt*k2/2), t+dt/2.)+self.sigma_r*self.__interpolate1D(eta,i+0.5)
				k33 = np.dot(jac,state+dt*k22/2)
				
				k4 = self.reaction((at+dt*k3), t+dt)+self.sigma_r*self.__interpolate1D(eta,i+1)
				k44 = np.dot(jac,state+dt*k33)
				
			at = at + dt*(k1+2*k2+2*k3+k4)/6.
			state = state + dt*(k11+2*k22+2*k33+k44)/6.
			#print(np.abs(at-last))
			deltas.append(np.abs(at-last)/dt)
			state, r = np.linalg.qr(state) 
			lyapSpec.append(np.abs(np.diag(r)))
			last = at
		return np.average(np.log(lyapSpec[-nit//4:]),axis=0)/dt, np.max(deltas)
	'''			
	
	def getNoisyLyapReaction(self,a0=None,beta=0,dt=0.1, nit=3000,eps=1e-5):

		lyap = []
		delta = 1e-6*(np.random.rand()-0.5)
		if a0 is None:
			at = self.a0+delta
		else:
			at = a0
		
		eta = cNoise.cNoise(beta=beta,shape=(nit+2,),std=1)
		eta = np.gradient(eta)
		
		lyapSpec=[]
		deltas = []
		for i in range(nit):
			t = i*dt
			angle = 2.0*np.pi*np.random.rand()
			at2 = at+eps*np.cos(angle)+eps*np.sin(angle)*1j
			if self.noiseType == 'multiplicative':
				k1 = self.reaction(at,t)+self.sigma_r*at*self.__interpolate1D(eta,i)
				k2 = self.reaction((at+dt*k1/2), t+dt/2.)+self.sigma_r*(at+dt*k1/2)*self.__interpolate1D(eta,i+0.5)
				k3 = self.reaction((at+dt*k2/2), t+dt/2.)+self.sigma_r*(at+dt*k2/2)*self.__interpolate1D(eta,i+0.5)
				k4 = self.reaction((at+dt*k3), t+dt)+self.sigma_r*(at+dt*k3)*self.__interpolate1D(eta,i+1)
			else:
				k1 = self.reaction(at,t)+self.sigma_r*self.__interpolate1D(eta,i)				
				k2 = self.reaction((at+dt*k1/2), t+dt/2.)+self.sigma_r*self.__interpolate1D(eta,i+0.5)
				k3 = self.reaction((at+dt*k2/2), t+dt/2.)+self.sigma_r*self.__interpolate1D(eta,i+0.5)
				k4 = self.reaction((at+dt*k3), t+dt)+self.sigma_r*self.__interpolate1D(eta,i+1)
			at = at + dt*(k1+2*k2+2*k3+k4)/6.
			
			if self.noiseType == 'multiplicative':
				k1 = self.reaction(at2,t)+self.sigma_r*at2*self.__interpolate1D(eta,i)
				k2 = self.reaction((at2+dt*k1/2), t+dt/2.)+self.sigma_r*(at2+dt*k1/2)*self.__interpolate1D(eta,i+0.5)
				k3 = self.reaction((at2+dt*k2/2), t+dt/2.)+self.sigma_r*(at2+dt*k2/2)*self.__interpolate1D(eta,i+0.5)
				k4 = self.reaction((at2+dt*k3), t+dt)+self.sigma_r*(at2+dt*k3)*self.__interpolate1D(eta,i+1)
			else:
				k1 = self.reaction(at2,t)+self.sigma_r*self.__interpolate1D(eta,i)				
				k2 = self.reaction((at2+dt*k1/2), t+dt/2.)+self.sigma_r*self.__interpolate1D(eta,i+0.5)
				k3 = self.reaction((at2+dt*k2/2), t+dt/2.)+self.sigma_r*self.__interpolate1D(eta,i+0.5)
				k4 = self.reaction((at2+dt*k3), t+dt)+self.sigma_r*self.__interpolate1D(eta,i+1)
			at2 = at2 + dt*(k1+2*k2+2*k3+k4)/6.
			
			deltas.append(np.abs(at-at2)/eps)
		return np.log(np.min(deltas))/dt, np.log(np.average(deltas))/dt, np.log(np.max(deltas))/dt		
	
	def __interpolate1D(self,noise,t):
		p1, p2 = int(np.floor(t)),int(np.ceil(t))
		if p1 == p2:
			return noise[p1]
		else:
			return noise[p1]*np.abs(t-p1)/np.abs(p1-p2) + noise[p2]*np.abs(t-p2)/np.abs(p1-p2)
	
	
	def getNoisyChainedSingleReaction(self,a0=None,beta=0,dt=0.1, nit=3000):
		'''
		Returns the iteration of a single amplitude (spatial part is ignored)
		
		The function integrates with rk4 method
		'''
		states = []
		delta = 1e-6*(np.random.rand()-0.5)
		if a0 is None:
			at = self.a0+delta
		else:
			at = a0
		
		eta = cNoise.cNoise(beta=beta,shape=(nit+2,),std=1)+1j*cNoise.cNoise(beta=beta,shape=(nit+2,),std=1)
		eta = np.gradient(eta)
		
		for i in range(nit):
			states.append(at)
			t = i*dt
			if self.noiseType == 'multiplicative':
				k1 = self.reaction(at,t)+self.sigma_r*at*self.__interpolate1D(eta,i)
				k2 = self.reaction((at+dt*k1/2), t+dt/2.)+self.sigma_r*(at+dt*k1/2)*self.__interpolate1D(eta,i+0.5)
				k3 = self.reaction((at+dt*k2/2), t+dt/2.)+self.sigma_r*(at+dt*k2/2)*self.__interpolate1D(eta,i+0.5)
				k4 = self.reaction((at+dt*k3), t+dt)+self.sigma_r*(at+dt*k3)*self.__interpolate1D(eta,i+1)
			else:
				k1 = self.reaction(at,t)+self.sigma_r*self.__interpolate1D(eta,i)
				k2 = self.reaction((at+dt*k1/2), t+dt/2.)+self.sigma_r*self.__interpolate1D(eta,i+0.5)
				k3 = self.reaction((at+dt*k2/2), t+dt/2.)+self.sigma_r*self.__interpolate1D(eta,i+0.5)
				k4 = self.reaction((at+dt*k3), t+dt)+self.sigma_r*self.__interpolate1D(eta,i+1)
			at = at + dt*(k1+2*k2+2*k3+k4)/6.
		return np.array(states)
		
	def reaction(self, a, t):
		a1 = a - (1+1j*self.c2)*(np.abs(a)**2)*a
		return np.array(a1)
		
	def interpolateNoise(self,time):
		'''
		Linear interpolation of the noise
		
		return the slice of nr and ni at the given time
		'''
		t = np.linspace(0,1,self.nr1.shape[0])
		p1, p2 = int(np.floor((self.nr1.shape[0]-1)*time)),int(np.ceil((self.nr1.shape[0]-1)*time))
		
		t1 = t[p1]
		t2 = t[p2]
		mr1 = self.nr1[p1]
		mr2 = self.nr1[p2]
		
		if np.abs(t1-t2)<1e-15:
			mr3 = mr1
		else:
			mr3 = np.abs(t1-time)*mr1/(np.abs(t2-t1))+np.abs(t2-time)*mr2/(np.abs(t2-t1))
		return mr3
		
		
	def solveRKF45(self,dt,ntimes,stepsave,dtTolerace=1e-4):
		state = self.getInitialCondition()
		times = []
		states = [state]	
			
		w = np.array([	[					0,0,0,0,0,0],
				[1/4,					0,0,0,0,0],
				[3/32,9/32,				0,0,0,0],
				[1932/2197,-7200/2197,7296/2197,	0,0,0],
				[439/216,-8,3680/513,-845/4104,	0,0],
				[-8/27, 2,-3544/2565,1859/4104,-11/40,0]
			])
		t = 0.0
		
		if 'beta' in self.noiseArgs:
			exponent = self.noiseArgs['beta']
		else:
			exponent = 2
		if 'std' in self.noiseArgs:
			std = self.noiseArgs['std']
		else:
			std = 0.01
		noiseShape = [ntimes]
		for i in range(self.dim):
			noiseShape.append(self.msize)
		self.nr1 = cNoise.cNoise(beta=exponent,shape=tuple(noiseShape),std=std)
		self.nr1 = np.gradient(self.nr1)[0]
			
		self.maxTime = (ntimes+2)*dt
				
		for time in tqdm.tqdm(range(ntimes)):
		
			step = dt
			
			k1 = step*self.timeDerivatives(state,								t		)
			k2 = step*self.timeDerivatives(state+k1*w[1,0], 		        				t+step/4	)
			k3 = step*self.timeDerivatives(state+k1*w[2,0]+k2*w[2,1], 					t+3*step/8	)
			k4 = step*self.timeDerivatives(state+k1*w[3,0]+k2*w[3,1]+k3*w[3,2],    				t+12*step/13	)
			k5 = step*self.timeDerivatives(state+k1*w[4,0]+k2*w[4,1]+k3*w[4,2]+k4*w[4,3],    		t+step		)
			k6 = step*self.timeDerivatives(state+k1*w[5,0]+k2*w[5,1]+k3*w[5,2]+k4*w[5,3]+k5*w[5,4],    	t+step/2	)
			
			approach4 = state + (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 -k5/5
			approach5 = state + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6
			
			error = np.max(np.abs(approach4-approach5))
			if error> dtTolerace:
				step = dt*((dtTolerace/(2*error))**.25)
			
				k1 = step*self.timeDerivatives(state,								t		)
				k2 = step*self.timeDerivatives(state+k1*w[1,0], 		        				t+step/4	)
				k3 = step*self.timeDerivatives(state+k1*w[2,0]+k2*w[2,1], 					t+3*step/8	)
				k4 = step*self.timeDerivatives(state+k1*w[3,0]+k2*w[3,1]+k3*w[3,2],    				t+12*step/13	)
				k5 = step*self.timeDerivatives(state+k1*w[4,0]+k2*w[4,1]+k3*w[4,2]+k4*w[4,3],    		t+step		)
				k6 = step*self.timeDerivatives(state+k1*w[5,0]+k2*w[5,1]+k3*w[5,2]+k4*w[5,3]+k5*w[5,4],    	t+step/2	)
				
				approach4 = state + (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 -k5/5
				
			t += step
			state = approach4 
			times.append(t)
			if time in stepsave:
				states.append(state)
		return np.array(states), np.array(times)
		
	def timeDerivatives(self,state,time):
		
		rFtState = fftn(np.real(state))
		iFtState = fftn(np.imag(state))
		normalizedTime = time/self.maxTime
		
		tnr1 = self.interpolateNoise(normalizedTime)
		
		# spectral shift variables
		specR = np.zeros(state.shape)
		specI = np.zeros(state.shape)
		
		# for every frequency dimension, updates the shift variable
		for i in range(len(state.shape)):
			
			# changing the dimension of the measured frequency since numpy's fftfreq is always 1D
			seq = np.ones(len(state.shape)).astype(int)
			fx  = 2*np.pi*fftfreq(state.shape[i])
			seq[i] = len(fx)
			fx = fx.reshape(*tuple(seq))
			
			#makes the shift
			specR = specR - (fx**2)*rFtState
			specI = specI - (fx**2)*iFtState
			
		# measures the laplacian from pseudospectral method
		lap  =    np.real(ifftn(specR)) + 1j*np.real(ifftn(specI))/(self.h**2)
		
		if self.noiseType == 'multiplicative':
			return (1+self.c1*1j)*np.array(lap) + self.reaction(state,time) + self.sigma_r*state*tnr1
		else:
			return (1+self.c1*1j)*np.array(lap) + self.reaction(state,time) + self.sigma_r*tnr1
