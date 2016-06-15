import numpy as np


class Model(object):
	def __init__(self,name,num):
		self.name=name
		self.num=num
	
	def make_mock(self,param):
		#np.random.seed()
                if self.name == "exp":
                        b = 1/param[0]
                        return np.random.exponential(b,self.num)
                if self.name == "normal":
                        sigm = np.diag(np.ones(len(param))*0.05) # how to generalise this?
                        return np.random.multivariate_normal(param,sigm,self.num)

