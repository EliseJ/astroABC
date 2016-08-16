import numpy as np


class Model(object):
	'''Toy class which simulates some data for testing astroABC sampler'''
	def __init__(self,name,num):
		self.name=name
		self.num=num
	
	def make_mock(self,param, var=0.05):
		'''
		Input: 
			param: variable to generate either exponential data or Normal data with variance = var
		'''
                if self.name == "exp":
                        b = 1/param[0]
                        return np.random.exponential(b,self.num)
                if self.name == "normal":
                        sigm = np.diag(np.ones(len(param))*var) 
                        return np.random.multivariate_normal(param,sigm,self.num)

