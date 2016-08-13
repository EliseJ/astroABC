import numpy as np
import os
import sys
sys.path.append(os.getcwd())
from .abc_class import *


def dist(d,x):
        '''Distance metric: rho'''
        return np.sum(np.abs(np.mean(x,axis=0) - np.mean(d,axis=0)))

class test_abc:

	def setUp(self):
		self.nparam =1
        	self.npart = 100
        	self.nsamples =1000
        	self.niter =15
        	self.tlevels = [.7,0.005]
        	self.model_type = "normal"

		self.prop={'tol_type':'exp',"verbose":1,'adapt_t':True,'threshold':75,
		'pert_kernel':2, 'variance_method':0, 'dist_type': "user",'dfunc':dist,
		'outfile':"mpi_test.txt",'mpi':False,'mp':False,'num_proc':None,
		'restart':"restart_test.txt",'from_restart':False}

		self.param = [0.1411]
		self.data = Model(self.model_type,self.nsamples).make_mock(self.param)

		priorname  = ["normal"]
		hyperp = [[x+0.2  ,0.05*2] for x in self.param]
		self.prior = zip(priorname,hyperp)


	def test_mpi(self):
		self.prop['mpi']=True
		sampler = ABC_class(self.nparam,self.npart,self.data,self.tlevels,self.niter,self.prior,**self.prop)
		assert(sampler.pool)
		


   



