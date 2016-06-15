import numpy as np


class Tolerance(object):
	def __init__(self,tol_type,tmin,tmax,nt):
		self.tol_type = tol_type
		self.nt=nt
		self.tmin = tmin ;self.tmax = tmax
		self.tol = self.set_tolerance()
	
	def set_tolerance(self):
		if self.tol_type =="const":
			return self.const_tol()
		elif self.tol_type =="linear":
			return self.linear_tol()
		elif self.tol_type =="exp":
			return self.exp_tol()
		elif self.tol_type =="log":
			return self.log_tol()
		else:
			print "Specify either const, linear, exp or log for tolerance class"

	def linear_tol(self):
		return np.linspace(self.tmax,self.tmin,num=self.nt)
		
	def log_tol(self):
		return np.logspace(self.tmax,self.tmin,num=self.nt)

	def const_tol(self):
		return np.ones(self.nt)*self.tmin

	def exp_tol(self):
		return np.logspace(np.log10(self.tmax), np.log10(self.tmin), num=self.nt)
