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
from setup_mpi_mp import *
import copy_reg
import types
from itertools import product

try:
        from mpi4py import MPI
        MPI = MPI
except ImportError:
        MPI = None

def wrapper_func(i):
	'''wrapper function to pass the correct sim pool to step() if mpi_splitcomm
	Input:
		tup_in: tuple (t = iteration number, Pid = particle id, weights, theta, variance, tol)
	Output:
		output from step() function
	'''
	if abcsampler.mpi_splitcomm:
		return step(i,abcsampler.parallel.sim_pool)
	else:
		return step(i)

def step(info_in,sim_pool=None):
                '''
                Function for a single particle to propose a new point in parameter space
                and accept it if rho(model(theta),data)<tol
                Input:
                        tup_in: tuple (t = iteration number, Pid = particle id)
                Returns:
                        trial_t - array of new parameters \theta
                        rho - distance rho(x,y)
                '''
                t, Pid,wgt,theta,variance,tol = info_in
                tm1=t-1

                while True:
                        if t ==0: #draw from prior
                                trial_t = [call_prior() for call_prior in abcsampler.prior]
                                if abcsampler.mpi_splitcomm:
                                        x = abcsampler.model(trial_t,sim_pool)
                                else:
                                        x = abcsampler.model(trial_t)
                                rho = abcsampler.dist(x)
                        else:
                                np.random.seed()
                                rpart = int(np.random.choice(abcsampler.npart,size=1,p=wgt[tm1]))
                                t_old = theta[tm1][rpart]
                                if abcsampler.variance_method == 4:
                                        covariance = variance[Pid]
                                else:
                                        covariance = variance
                                trial_t = np.atleast_1d(scipy.stats.multivariate_normal.rvs(mean= t_old,cov=covariance,size=1))
                                if abcsampler.mpi_splitcomm:
                                        x = abcsampler.model(trial_t,sim_pool)
                                else:
                                        x = abcsampler.model(trial_t)
                                rho = abcsampler.dist(x)
                        if rho <= tol[t]:
                                break

                return trial_t,rho


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
			tol_type: tolerance level setting. Can be "linear","log", "exp", "const". 
			verbose: 0 = no output to screen, 1  = print out to screen
			adapt_t: Boolean True/False for adaptive threshold setting
			threshold: quantile level if adapt_t == True.
			pert_kernel: 1 =  component wise perturbation with local diag variance; 2 = multivariate perturbation based on local covariance
			variance_method: 0=Weighted covariance, 1=Filippi et al 2012, eq 12 & 13, 2=Simple variance estimate (Turner & Van Zandt 2012), 3=Leodoit_Wolf, 4=k-nn
			k_near: int, number of nearest neighbours for local covariance
			dist_type: string, distance metric method. Default is "user"
			dfunc: method name when dist_type == user"
			datacov: string, data covariance file if needed in distance metric
			outfile: string, name of output file 
			mpi: Boolean, True/False
			mpi_splitcomm: Boolean, True/False
			num_abc: int, number of procsessors for abc particles if split_mpicomm == True
			mp: Boolean, True/False
			num_proc: int, if mp=True
			restart: string, name of restart file
			from_restart: Boolean, True/False
		'''
		prop_defaults={"tol_type":"exp","verbose":0,'adapt_t':False,
		'threshold':75,'pert_kernel':1,'variance_method':0,'k_near':5,'dist_type': "user",         
		'dfunc':None,'datacov':None,'outfile':'abc_out.txt','mpi': None, 'mp':None,'num_proc':None,'mpi_splitcomm':False,
		'num_abc':None,'restart':None,'from_restart':False}
		for (prop, default) in prop_defaults.iteritems():
			setattr(self, prop, kwargs.get(prop, default))

		if self.from_restart:
                	backup_files(self.outfile)

		global abcsampler
        	abcsampler = self  

		#setup mpi or mp
		self.parallel = Parallel(self.mpi, self.mp,self.mpi_splitcomm, self.num_proc,self.num_abc,self.verbose)

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
		elif self.variance_method == 4:
			self.Variance = nearest_neighbour(nparam,npart,self.pert_kernel)         
		else: #default
			self.Variance = weighted_cov(nparam,npart,self.pert_kernel) 

		self.unpack_priors(priors)
		self.end_sampling = False
		if self.verbose and self.parallel.master==0:
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


	def sample(self,model_simulator):
		'''
		Begin sampling  
		Input:
		 	model_simulator  - func which simulates data
		'''
		self.model=model_simulator
		
		#if self.parallel.rank ==0:print "\t Running sampler...\t "

		if self.from_restart:
			t,th,wgt,dist=read_restart_files(self.restart,self.nparam,self.npart)
			self.theta[t] = th ; self.wgt[t]=wgt; self.Delta[t] = dist
			if self.adapt_t and t <self.niter-1:
                         	self.tol[t+1]=self.iteratively_adapt(t)
			ctr = t+1
		else:
			ctr=0

		while self.end_sampling == False:
			if self.mpi_splitcomm: 
				if self.parallel.rank in self.parallel.abc_ranks:
					ctr = self.sample_loop(ctr)
				else: #worker node for sim, wait until sim pool is closed 
					self.parallel.sim_pool.worker()
					self.end_sampling = True
			else:
				ctr = self.sample_loop(ctr)

		if self.mpi or self.mp:
			if self.mpi_splitcomm:
				if self.parallel.rank in self.parallel.abc_ranks:
					self.parallel.pool.close()
					if self.parallel.rank !=0:
						self.parallel.sim_pool.close()
			else:
				self.parallel.pool.close()
		
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

			self.variance = self.calculate_variance(t)
		else:
			self.variance =0

		if self.mpi or self.mp:
			if self.mp:
				pool_outputs = self.parallel.pool.map(self.classstep, zip([t]*(self.npart),range(self.npart)))
			else:
				pool_outputs = self.parallel.pool.map(wrapper_func, zip([t]*(self.npart),\
				range(self.npart),[self.wgt]*(self.npart),\
				[self.theta]*(self.npart),[self.variance]*(self.npart), [self.tol]*(self.npart)))

			for i in range(self.npart):
				if pool_outputs: # prevent error when mpi worker pool is closed 
					self.theta[t][i],self.Delta[t][i]  = pool_outputs[i] 
				else:
					self.end_sampling = True
					return  #pool is closed so worker just returns
			if t:
				self.wgt[t] = self.parallel.pool.map(self.particle_weight, zip([t]*(self.npart),range(self.npart)))
		else:
			for i in np.arange(self.npart):                                          
				self.theta[t][i],self.Delta[t][i] = step((t,i, self.wgt, self.theta, self.variance, self.tol))			
				if t:
					self.wgt[t][i] = self.particle_weight((t,i))
		#normalize
		self.wgt[t] = self.wgt[t]/np.sum(self.wgt[t])

			
                if self.outfile and self.parallel.master==0:
                	write_to_file(t,self.outfile,self.nparam,self.npart,self.theta[t],self.Delta[t],self.wgt[t])
			if self.restart and t:
				write_restart_file(self.restart,t,self.theta[t],self.wgt[t],self.Delta[t],self.nparam,self.npart)
                sys.stdout.flush()
		
		if self.verbose and self.parallel.master==0:
			print "\t Step:",t,"\t tol:",self.tol[t],"\t Params:",[np.mean(self.theta[t][:,ii]) for ii in range(self.nparam)]
		if self.adapt_t and t <self.niter-1:
			self.tol[t+1]=self.iteratively_adapt(t)
		return t+1



	def calculate_variance(self,t):
		if self.variance_method ==1:
			return self.Variance.get_var(t,self.theta[t-1],self.Delta[t-1],self.tol[t-1],self.wgt[t-1])
		elif self.variance_method ==2 or self.variance_method ==3:
			return self.Variance.get_var(t,self.theta[t-1])
		elif self.variance_method == 4:
			return self.Variance.get_var(t,self.theta[t-1],self.wgt[t-1],self.k_near)
		else:
			return self.Variance.get_var(t,self.theta[t-1],self.wgt[t-1])
		
	def iteratively_adapt(self,t):
		'''
		Drovandi & Pettitt 2011, use qth quantileof the distance for t iterations
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
		if self.variance_method == 4:
		    covariance = self.variance[Pid]
		else:
		    covariance = self.variance
                if np.linalg.det(covariance) <1.E-15:
                    #maybe singular matrix; check diagonals for small values
		    #if self.verbose:
                    #	print "Variance is a singular matrix", covariance
                    #	print "using l2 shrinkage with the Ledoit-Wolf estimator..."
                    covariance, _  =  ledoit_wolf(self.theta[t])
                return scipy.stats.multivariate_normal(mean=self.theta[t][Pid],cov=covariance,allow_singular=True).pdf

	def classstep(self,info_in):
                '''
		if mp==True.
		copy of step function above; messy but necessary for mp which can't access global abcsampler -EJ
                '''
                t, Pid = info_in
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
                                if self.variance_method == 4:
                                        covariance = self.variance[Pid]
                                else:
                                        covariance = self.variance
                                trial_t = np.atleast_1d(scipy.stats.multivariate_normal.rvs(mean= t_old,cov=covariance,size=1))
                                x = self.model(trial_t)
                                rho = self.dist(x)
                        if rho <= self.tol[t]:
                                break
                return trial_t,rho
