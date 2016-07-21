import numpy as np
import subprocess

def write_to_file(t,outfile,nparam,npart,theta,delta,wgt):
	'''Method to write sampler outputs to txt file
	Input: 
	t - iteration level
	outfile - str, filename
	nparam - int, number of parameters
	npart - int, number of particles
	theta - arr, sampler outputs at iteration t
	delta -arr, distances at iteration t
	wgt - arr, weights at iteration t
	'''
        if t==0:
                f = open(outfile,'w',0)
                for ncount in range(nparam):
                        f.write("param#%s \t "%ncount)
                f.write("dist \t  wgt \n")
        else:
                f = open(outfile,'a')

        data= theta.flatten().reshape(npart,nparam)
        dists= delta.flatten().reshape(npart,1)
        wgts= wgt.flatten().reshape(npart,1)
        output = np.hstack((data,dists,wgts))
        for o in output:
                for indv in o:
                        f.write("%f \t" % indv)
                f.write("\n")
        f.flush()
        f.close()

def write_restart_file(restart,t,theta,wgt,dist,nparam,npart):
	'''Method to write restart txt file
        Input: 
        restart - str, filename
        t - iteration level
        nparam - int, number of parameters
        npart - int, number of particles
        theta - arr, sampler outputs at iteration t
        delta -arr, distances at iteration t
        wgt - arr, weights at iteration t
        '''

	if t>1:backup_files(restart)
	data= theta.flatten().reshape(npart,nparam)
        wgts= wgt.flatten().reshape(npart,1)
        dists= dist.flatten().reshape(npart,1)
        output = np.hstack((data,dists,wgts))
	f = open(restart,"w",0)
	f.write("%f\n"%t)
	for o in output:
                for indv in o:
                        f.write("%.8f \t" % indv)
                f.write("\n")
        f.flush()
        f.close()

def backup_files(filename):
	'''Method to backup file before writing any output
        Input: 
        filename - str, filename
	'''
	cmd = 'cp  '+ filename+'  '+filename+'.bck'
	try:
		process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
	except:
		pass


def read_restart_files(restart,nparam,npart):	
	'''Method to read restart txt file
        Input: 
        restart - str, filename
        nparam - int, number of parameters
        npart - int, number of particles
	Returns:
        t - iteration level
        theta - arr, sampler outputs at iteration t
        wgt - arr, weights at iteration t
        dist -arr, distances at iteration t
        '''
	print "\t -----> Reading restart files .... \t \n"
	t = np.genfromtxt(restart,max_rows=1)
	data=np.genfromtxt(restart,skip_header=1)
	theta = data[:,:-2].reshape(npart,nparam)
	wgt = data[:,-1]
	dist = data[:,-2]
	#make sure normalised weights
	wgt = wgt/np.sum(wgt)
	return int(t),theta,wgt,dist
		
