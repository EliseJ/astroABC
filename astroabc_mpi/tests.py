import numpy as np
from .abc_class import *


def dist(d,x):
        '''Distance metric: rho'''
        return np.sum(np.abs(np.mean(x,axis=0) - np.mean(d,axis=0)))

class test_abc:

	def setUp(self):
		self.nparam =2
        	self.npart = 100
        	self.nsamples =1000
        	self.niter =15
        	self.tlevels = [.7,0.005]
        	self.model_type = "normal"

		self.prop={'tol_type':'exp',"verbose":1,'adapt_t':True,'threshold':75,
		'pert_kernel':2, 'variance_method':0, 'dist_type': "user",'dfunc':dist,
		'outfile':"mpi_test.txt",'mpi':False,'mp':False,'num_proc':None,
		'restart':"restart_test.txt",'from_restart':False}

		self.param = [0.1411,0.8723]
		self.data = Model(self.model_type,self.nsamples).make_mock(self.param)

		priorname  = ["normal", "normal"]
		hyperp = [[x+0.2  ,0.05*2] for x in self.param]
		self.prior = zip(priorname,hyperp)


	def test(self):
		sampler = ABC_class(self.nparam,self.npart,self.data,self.tlevels,self.niter,self.prior,**self.prop)
		for i in range(self.niter):
			param_means = [np.mean(sampler.theta[i][:,j]) for j in range(self.nparam)]
			for p in range(self.nparam):
				assert(param_means[p] < np.inf)


	def test_tolerance(self):
		tol = Tolerance('exp',self.tlevels[1],self.tlevels[0],self.niter).tol
		for i in range(len(tol)-1):
			assert(tol[i] > tol[i+1])

		tol = Tolerance('const',self.tlevels[1],self.tlevels[0],self.niter).tol
		for i in range(len(tol)-1):
			assert(tol[i] == tol[i+1])

		tol = Tolerance('linear',self.tlevels[1],self.tlevels[0],self.niter).tol
		for i in range(len(tol)-1):
			assert(tol[i] > tol[i+1])
		
		tol = Tolerance('log',self.tlevels[1],self.tlevels[0],self.niter).tol
		for i in range(len(tol)-1):
			assert(tol[i] > tol[i+1])

		
	def test_mp(self):
		self.prop['mp']=True
		self.prop['num_proc']=2
		sampler = ABC_class(self.nparam,self.npart,self.data,self.tlevels,self.niter,self.prior,**self.prop)
		assert(sampler.pool)

	def test_variance(self):
		sampler = ABC_class(self.nparam,self.npart,self.data,self.tlevels,1,self.prior,**self.prop)
		for p1 in self.nparam:
			for p2 in self.nparam:
				assert(sampler.variance[p1][p2] < np.inf)
		theta = sampler.theta
		
		


   



