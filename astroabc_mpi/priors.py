import numpy as np
import scipy.stats
import time


class Prior_class(object):
        '''priors to choose from , prior is array of keywords ["uniform, "gauss"] for each param"'''
        def __init__(self,priorname,hyperparams):
                self.priorname=priorname
                self.hyperparams = hyperparams
        
        def return_priorprob(self,value):
		#np.random.seed()
                if self.priorname =="gamma":
                        x = 1./self.hyperparams[1]
                        return scipy.stats.gamma.pdf(value, self.hyperparams[0],scale=x)
		elif self.priorname =="normal":
                        return scipy.stats.norm.pdf(value, loc = self.hyperparams[0],scale=self.hyperparams[1])
		elif self.priorname =="uniform":
			width = self.hyperparams[1] - self.hyperparams[0]
                        return scipy.stats.uniform.pdf(value, loc = self.hyperparams[0],scale=width)


        def prior(self):
		np.random.seed()
                if self.priorname =="gamma":
                        k=self.hyperparams[0]
                        scale = 1./self.hyperparams[1]
                        return float(np.random.gamma(k,scale))
                elif self.priorname =="normal":
                        return float(np.random.normal(self.hyperparams[0],self.hyperparams[1],size=1))
                elif self.priorname =="uniform":
                        return float(np.random.uniform(low=self.hyperparams[0],high=self.hyperparams[1],size=1))

