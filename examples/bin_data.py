import numpy as np
import math

def read_fitres(fname):
                dtype=[('z',np.float32),('mu',np.float32),('mu_err',np.float32)]
                data=np.loadtxt(fname,skiprows = 15,dtype =dtype, usecols=(34,48,50))
                return data

def bin_data(z,mu,err,Nbins=10):
    binw = (z.max()-z.min())/Nbins
    #print binw
    zbins = np.arange(z.min(),z.max(),binw)
    mu_in_bin=[];err_in_bin=[]
    avmu_bin=[] ; averr_bin=[]
    for i in range(Nbins):
        mu_in_bin.append([]); err_in_bin.append([])
    for i,m in enumerate(mu):
        for j in range(Nbins): 
            if z[i] >=z.min()+ j*binw and z[i] < z.min() + (j+1)*binw:
                mu_in_bin[j].append(m) ; err_in_bin[j].append(err[i])
    for i in range(Nbins):
        if mu_in_bin[i]:
            avmu_bin.append(np.mean(mu_in_bin[i]))
        else: avmu_bin.append(0)
        if err_in_bin[i]:
            averr_bin.append(np.mean(err_in_bin[i]))
        else: averr_bin.append(0)
    return  zbins, avmu_bin,averr_bin,mu_in_bin,mu_in_bin

def read_data():
	data = read_fitres("example_data_file.txt")

	#Keep data points with z>0.5
	z_new = data['z'][0]
	mu_new = data['mu'][0]
	err_new = data['mu_err'][0]

	for i in range(1,345):
	    if data['z'][i] >= 0.5:
        	z_new = np.append(z_new,data['z'][i])
        	mu_new = np.append(mu_new,data['mu'][i])
		err_new = np.append(err_new,data['mu_err'][i])

	#bin this data
	zbins,avmu_bin,averr_bin,mu_in_bin_new,mu_in_bin_new = bin_data(z_new,mu_new,err_new,Nbins=20)
	return zbins,avmu_bin,averr_bin,mu_in_bin_new,mu_in_bin_new


