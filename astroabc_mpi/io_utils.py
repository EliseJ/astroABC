import numpy as np
import subprocess

def write_to_file(t,outfile,nparam,npart,theta,delta,wgt):
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
	cmd = 'cp  '+ filename+'  '+filename+'.bck'
	try:
		process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
	except:
		pass


def read_restart_files(restart,nparam,npart):	
	print "\t -----> Reading restart files .... \t \n"
	t = np.genfromtxt(restart,max_rows=1)
	data=np.genfromtxt(restart,skip_header=1)
	theta = data[:,:-2].reshape(npart,nparam)
	wgt = data[:,-1]
	dist = data[:,-2]
	#make sure normalised weights
	wgt = wgt/np.sum(wgt)
	return int(t),theta,wgt,dist
		
