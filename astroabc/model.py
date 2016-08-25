import numpy as np


class Model(object):
	'''Toy class which simulates some data for testing astroABC sampler'''
	def __init__(self,name,num):
		self.name=name
		self.num=num
	
	def make_mock(self,(param, var)):
		'''
		Input: 
			param: variable to generate either exponential data or Normal data with variance = var
		'''
                if self.name == "exp":
                        b = 1/param[0]
                        return np.random.exponential(b,self.num)
                if self.name == "normal":
                        if isinstance(var,float):
				sigm = np.diag(np.ones(len(param))*var) 
			elif len(var)==len(param):
				sigm = np.diag(var)
			else:
				sigm = var.reshape(len(param),len(param)) #covariance matrix
                        return np.random.multivariate_normal(param,sigm,self.num)

