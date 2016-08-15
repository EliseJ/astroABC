import numpy as np
import scipy.stats
import time

NFILES = 5
FILENAME = "omw_" # prefix for files followed by filenumber

class Prior_class(object):
        '''Prior class''' 
	
        def __init__(self,priorname,hyperparams):
		'''Input:
		priorname - array of keywords ["uniform, "gauss"] for each param
		hyperparams - array of arrays [[min,max], [mu, variance],...] for each param
		'''
                self.priorname=priorname
                self.hyperparams = hyperparams
		if self.priorname == "nonstandard": #first hyperparam = column from file to be read, specify filename and number of files above
			'''useful for sampling from non-standard discrete pdf e.g. Planck/WMAP chain'''
			self.read_data(hyperparams[0])
			self.inv_transform_spline()
			self.pdf_spline()

	def read_data(self,colnum):
		'''Only for "nonstandard". Method to read discrete pdf for parameter from file
		Input: colnum: column number to be read from file
		'''
		self.param=[]
		for i in range(1,NFILES):
			d = np.loadtxt(FILENAME+str(i)+".txt")
			for j in range(len(d[:,colnum])):
				self.param.append(d[:,colnum][j])

	def inv_transform_spline(self):
		'''Only for "nonstandard". Method to create inverse spline to discrete cumulative distribution function 
		to allow drawing random variables.
		Warning: user should check that spline faithfully matches actual cdf.
		'''
		srt_param=np.sort(self.param)
		cdf = np.array(range(len(self.param)))/float(len(self.param))
		#create a spline
		self.spline2_cdf = UnivariateSpline(cdf,srt_param,k=5)

	def pdf_spline(self):
		'''Only for "nonstandard". Method creates a spline to the normalised PDF for discrete parameter values.
		Warning: user should check that spline faithfully matches actual pdf.
		'''
		hist,nbins = np.histogram(self.param,normed=True,bins=200)
		self.spline2_pdf = interp1d(nbins[1:],hist)
        
        def return_priorprob(self,value):
		'''Input:
		value -  random variable
		Returns:
		probability of rv given the prior dist
		'''
                if self.priorname =="gamma":
                        x = 1./self.hyperparams[1]
                        return scipy.stats.gamma.pdf(value, self.hyperparams[0],scale=x)
		elif self.priorname =="normal":
                        return scipy.stats.norm.pdf(value, loc = self.hyperparams[0],scale=self.hyperparams[1])
		elif self.priorname =="uniform":
			width = self.hyperparams[1] - self.hyperparams[0]
                        return scipy.stats.uniform.pdf(value, loc = self.hyperparams[0],scale=width)
		elif self.priorname == "nonstandard":
			return self.spline2_pdf(value)


        def prior(self):
		'''
		Returns a random variable from the prior distribution
		'''
		np.random.seed()
                if self.priorname =="gamma":
                        k=self.hyperparams[0]
                        scale = 1./self.hyperparams[1]
                        return float(np.random.gamma(k,scale))
                elif self.priorname =="normal":
                        return float(np.random.normal(self.hyperparams[0],self.hyperparams[1],size=1))
                elif self.priorname =="uniform":
                        return float(np.random.uniform(low=self.hyperparams[0],high=self.hyperparams[1],size=1))
		elif self.priorname == "nonstandard":
			uni_rvs = np.random.uniform()
			return float(self.spline2_cdf(uni_rvs))



