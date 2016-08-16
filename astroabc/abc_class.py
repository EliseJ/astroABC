############### astro ABC #######
#####author: Elise Jennings #####
####elise@fnal.gov ##########

import numpy as np
import scipy.stats
import sys
from sklearn.covariance import GraphLassoCV, ledoit_wolf
import os
sys.path.append(os.getcwd())

from myutils import *
from tolerance import *
from variance import *
from priors import *
from model import *
from io_utils import *
from mpi_pool import MpiPool
import multiprocessing as mp
import copy_reg
import types
from itertools import product



class ABC_class(object):
	'''Approximate Bayesian Computation Sequentual Monte Carlo class'''
	def __init__(self,nparam,npart,data,tlevels,niter,priors,**kwargs):
		'''
		Input:
			nparam: number of parameters to vary in sampler
			npart: number of particles/walkers for Sequential Monte Carlo sampler
			data: input data vector
			tlevels: [max,min] tolerance levels max and min values
			niter: number of iteration levels
			priors: list of tuples (priorname, [hyperparameters for prior])
			kwargs: dictionary of key words ; defaults given below
		'''
		prop_defaults={"tol_type":"exp","verbose":0,'adapt_t':False,
		'threshold':75,'pert_kernel':1,'variance_method':0, 'dist_type': "user",
		'dfunc':None,'datacov':None,'outfile':'abc_out.txt','mpi': None, 'mp':None,'num_proc':None,
		'restart':None,'from_restart':False}
		for (prop, default) in prop_defaults.iteritems():
			setattr(self, prop, kwargs.get(prop, default))


		if self.from_restart:
                	backup_files(self.outfile)
                if self.mpi:
                        self.pool = MpiPool()
			self.master = self.pool.rank
			if self.pool.size ==1:
				print "\n \t Please run using mpirun -np #proc if you want to run using mpi. \n\texiting ..."
				sys.exit(0)
		elif self.mp:
			self.pool = mp.Pool(self.num_proc)
			self.master  = 0
		else:
			self.master =0
		check_input(nparam,npart,priors,tlevels[0],tlevels[1],self.dist_type,self.datacov,self.dfunc)
		
		self.data = data
		self.npart = npart
		self.nparam = nparam
		self.niter = niter
		self.tmin = tlevels[1]

		self.allocate()
		self.tol = Tolerance(self.tol_type,tlevels[1],tlevels[0],niter).tol

		if self.variance_method ==1:
			self.Variance = Filippi(nparam,npart,self.pert_kernel) 
		elif self.variance_method ==2:
			self.Variance = TVZ(nparam,self.pert_kernel) 
		elif self.variance_method ==3:
			self.Variance = Leodoit_Wolf(nparam,self.pert_kernel) 
		else: #default
			self.Variance = weighted_cov(nparam,npart,self.pert_kernel) 

		self.unpack_priors(priors)
		self.end_sampling = False
		if self.verbose and self.master==0:
			print_header(npart,niter,self.tol_type,tlevels,priors)
		

	def unpack_priors(self,priors):
		'''
		Parameters of interest, theta, may have different priors from Priors class
		input: 
			priors = list of ['name of prior', [hyperparams of prior]] for each t in theta, 
			currently, gamma, uniform and normal supported.
		returns: 
			methods from instances of Prior class which allow sample to generate rvs 
			using self.prior() and  pdf(number) using self.priorprob(number)
		'''
		self.prior = np.empty(self.nparam,dtype=np.object)
		self.priorprob = np.empty(self.nparam,dtype=np.object)
		for i,p in enumerate(priors):
			pcls = Prior_class(p[0],p[1])
			self.prior[i] = np.vectorize(pcls.prior)
			self.priorprob[i]=np.vectorize(pcls.return_priorprob)
	def allocate(self):
		'''allocate numpy arrays for parameter values, weights and distances'''
		self.theta=np.zeros([self.niter,self.npart,self.nparam])	
		self.wgt=np.zeros([self.niter,self.npart])
		self.Delta=np.zeros([self.niter,self.npart])

	def dist(self,x):
		'''distance metric function
		Input:
		x - simulationed data sample at proposed parameter value
		Returns:
		distance to be compared to the threshold at iteration t
		'''
        	if self.dist_type == "mean":
			return np.sum(np.abs(np.mean(x,axis=0) - np.mean(self.data,axis=0)))
        	if self.dist_type == "chi2":
			d=x-self.data
			return np.dot(d,np.dot(self.datacov,d))
		if self.dist_type == "user":
			return self.dfunc(self.data,x)



	def step(self,t):
		'''
		Method for a single particle to propose a new point in parameter space
		and accept it if rho(model(theta),data)<tol
		Input:
			t - iteration number
		Returns:
			trial_t - array of new parameters \theta
			rho - distance rho(x,y)
		'''
		tm1=t-1
		while True:
			if t ==0: #draw from prior
				trial_t = [call_prior() for call_prior in self.prior]
				x = self.model(trial_t)
				rho = self.dist(x)
			else:
                                np.random.seed()
				rpart = int(np.random.choice(self.npart,size=1,p=self.wgt[tm1]))
				t_old = self.theta[tm1][rpart]		
				trial_t = np.atleast_1d(scipy.stats.multivariate_normal.rvs(mean= t_old,cov=self.variance,size=1))
				x = self.model(trial_t)
				rho = self.dist(x)
			if rho <= self.tol[t]:
				break
		return trial_t,rho


	def sample(self,model_simulator):
		'''
		Begin sampling  
		Input:
		 	model_simulator  - func which simulates data
		'''
		self.model=model_simulator
		
		print "\t Running sampler...\t "

		if self.from_restart:
			t,th,wgt,dist=read_restart_files(self.restart,self.nparam,self.npart)
			self.theta[t] = th ; self.wgt[t]=wgt; self.Delta[t] = dist
			if self.adapt_t and t <self.niter-1:
                         	self.tol[t+1]=self.iteratively_adapt(t)
			ctr = t+1
		else:
			ctr=0

		while self.end_sampling == False:
			ctr = self.sample_loop(ctr)
		if self.mpi or self.mp:
			self.pool.close()
		
	def sample_loop(self,t):
		'''
		At each iteration t: 
			-each particle finds a point in parameter space that satifies rho(model(theta),data)<tol
			-weights are calculated for each particle
			-variances of parameters calculated
			-weights are normalized
		input:
			iter t
		'''
		if  t+1 == self.niter or self.tol[t] == self.tmin:
			self.end_sampling = True

		if not(t): 
			self.wgt[t] =1./self.npart
		
		if t: 	
			if self.variance_method ==1:
				self.variance = self.Variance.get_var(t,self.theta[t-1],self.Delta[t-1],self.tol[t-1],self.wgt[t-1])
			elif self.variance_method ==2 or self.variance_method ==3:
				self.variance = self.Variance.get_var(t,self.theta[t-1])
			else:
				self.variance = self.Variance.get_var(t,self.theta[t-1],self.wgt[t-1])
		if self.mpi or self.mp:
        		pool_outputs = self.pool.map(self.step, [t]*(self.npart))
			for i in range(self.npart):
				if pool_outputs: # prevent error when mpi worker pool is closed 
					self.theta[t][i],self.Delta[t][i]  = pool_outputs[i] 
				else:
					self.end_sampling = True
					return  #pool is closed so worker just returns
			if t:
				self.wgt[t] = self.pool.map(self.particle_weight, zip([t]*(self.npart),range(self.npart)))
		else:
			for i in np.arange(self.npart):
				self.theta[t][i],self.Delta[t][i] = self.step(t)
				if t:
					self.wgt[t][i] = self.particle_weight((t,i))
		#normalize
		self.wgt[t] = self.wgt[t]/np.sum(self.wgt[t])

			
                if self.outfile and self.master==0:
                	write_to_file(t,self.outfile,self.nparam,self.npart,self.theta[t],self.Delta[t],self.wgt[t])
			if self.restart and t:
				write_restart_file(self.restart,t,self.theta[t],self.wgt[t],self.Delta[t],self.nparam,self.npart)
                sys.stdout.flush()
		
		if self.verbose and self.master==0:
			print "\t Step:",t,"\t tol:",self.tol[t],"\t Params:",[np.mean(self.theta[t][:,ii]) for ii in range(self.nparam)]
		if self.adapt_t and t <self.niter-1:
			self.tol[t+1]=self.iteratively_adapt(t)
		return t+1

		
	def iteratively_adapt(self,t):
		'''
		Turner & Van Zandt (2012), use qth quantileof the distance for t iterations
		unless we hit the tmin requested
		'''
		new_tol= np.percentile(self.Delta[t], self.threshold)
		if new_tol < self.tmin: 
			new_tol = self.tmin
		return new_tol


	def particle_weight(self,tup_in):
		'''
		At each iteraction, t, this method calculates the weight of particle i
                at t not equal to 0 weight is calculated according to kernel. 
                input: 
                        tup_in is a tuple of (iter t,particle id Pid)
                '''
		t, Pid = tup_in
                tm1 = t-1

                # apply kernel to each particle's parameter vector
                Kf = self.kernel(Pid,t)
                kernels = Kf(self.theta[t-1])
                if  np.any(self.wgt[tm1]) ==0 or np.any(kernels)==0:
                        print "Error computing Kernel or weights...", kernels, self.wgt[tm1]
                        sys.exit(1)
                priorproduct = np.prod([f[0](f[1]) for f in np.vstack((self.priorprob,self.theta[t][Pid])).T])
		return priorproduct/(np.sum(self.wgt[tm1]*kernels))


        def kernel(self,Pid,t):
                if np.linalg.det(self.variance) <1.E-15:
                    #maybe singular matrix; check diagonals for small values
                    print "Variance is a singular matrix", self.variance
                    print "using l2 shrinkage with the Ledoit-Wolf estimator..."
                    self.variance, _  =  ledoit_wolf(self.theta[t])
                return scipy.stats.multivariate_normal(mean=self.theta[t][Pid],cov=self.variance).pdf

	#http://stackoverflow.com/questions/25382455/python-notimplementederror-pool-objects-cannot-be-passed-between-processes
	def __getstate__(self):
		self_dict = self.__dict__.copy()
		del self_dict['pool']
		return self_dict

	def __setstate__(self, state):
		self.__dict__.update(state)

	
def _pickle_method(method):
        """
        http://stackoverflow.com/questions/11726809/
        python-efficient-workaround-for-multiprocessing-a-function-that-is-a-data-member
        """
        func_name = method.im_func.__name__
        obj = method.im_self
        cls = method.im_class
        cls_name = ''
        if func_name.startswith('__') and not func_name.endswith('__'):
                cls_name = cls.__name__.lstrip('_')
        if cls_name:
                func_name = '_' + cls_name + func_name
        return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
        for cls in cls.mro():
                try:
                        func = cls.__dict__[func_name]
                except KeyError:
                        pass
                else:
                        break
        return func.__get__(obj, cls)

copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
