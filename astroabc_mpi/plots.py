from itertools import combinations
import matplotlib.pylab as plt
import seaborn as sns
import sys

class Plotter(object):
	def __init__(self,sampler,true,plot_iter=False):
		self.splr=sampler
		self.t = sampler.theta
		self.true=true
		self.plot_iter=plot_iter

	def plot(self):
        	if self.splr.nparam==1:
                	#true_p = self.splr.true_post(500,0.1,0.1)
                	ax=sns.distplot(self.t.flatten(), color='green', hist=True,label="ABC Estimated posterior")
                	ax.plot([self.true[0],self.true[0]],[0,10],color='red',label="True")
                	#ax.hist(self.t.flatten(),color="cyan",normed=1, histtype='bar',label=" ABC Estimated posterior")
                	ax.legend(frameon=False,loc=2)
                	ax.set_xlabel('l')
                	plt.savefig("example_plot.png")
        	else:
                	for p in combinations(range(self.splr.nparam),2):
                        	x=self.t.flatten().reshape(self.splr.npart*self.splr.niter,self.splr.nparam).T[p[0]]
                        	y=self.t.flatten().reshape(self.splr.npart*self.splr.niter,self.splr.nparam).T[p[1]]
				if self.plot_iter:
					for i in range(self.splr.niter):
						xt= x[self.splr.npart*i:(self.splr.npart-1)*(i+1)]
						yt=y[self.splr.npart*i:(self.splr.npart-1)*(i+1)]
						try:
							self.plot_chain(xt,yt,p,i)
						except: print("Error plotting KDE contour!"); sys.exit(0)
					try:
						self.plot_chain(x,y,p,"_tot")
					except: print("Error plotting KDE contour!"); sys.exit(0)
				else:
					try:
						self.plot_chain(x,y,p,"_tot")
					except: print("Error plotting KDE contour!"); sys.exit(0)


	def plot_chain(self,x,y,p,i):
		sns.plt.xlim(-1.,1.)
		sns.plt.ylim(-1.,1.)
		ax = sns.kdeplot(x, y, shade=True)
		plt.plot(self.true[p[0]],self.true[p[1]],"o",color="red",markersize=3)
		plt.savefig("./plots/example_plot"+str(p[0])+str(p[1])+"_iter"+str(i)+".png")
		plt.close()	
