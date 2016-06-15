import numpy as np
import scipy.stats

from abc_class import ABC_class

class Weights(ABC_class):
	def __init__(self):
		pass

	def particle_weight(self,t,i):
                '''
                At each iteraction, t, this method calculates the weight of particle i
                at t=0 returns w = 1/npart for each particle
                at t not equal to 0 weight is calculated according to kernal. 
                input: 
                        iter t
                        particle id i
                '''
                tm1 = t-1

                # apply kernal to each particle's parameter vector
                Kf = self.kernal(i,t)
                kernals = Kf(self.theta[tm1])

                if  np.any(self.wgt) ==0 or np.any(kernals)==0:
                        print "Error computing Kernal or weights..."
                        sys.exit(1)
                priorproduct = np.prod([f[0](f[1]) for f in np.vstack((self.priorprob,self.theta[t][i])).T])
                self.wgt[t][i] = priorproduct/(np.sum(self.wgt[tm1]*kernals))


        def kernal(self,i,t):
                #component wise pertubation
                cov = np.diag(self.var[t]) # this is the part where you might just use KNN to calculate kernal with weighted cov
                return scipy.stats.multivariate_normal(mean=self.theta[t][i],cov=cov).pdf
