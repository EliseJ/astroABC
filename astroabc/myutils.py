import sys

def check_input(nparam,npart,priors,t1,t2,dtype,datacov,dfunc):
	'''Simple error checking on input before we start sampling
	Input:
	nparam - number of paramerer
	npart - number of particles
	priors - priots
	t1 - max treshold
	t2 - min threshold
	dtype - data covariance matrix, this toy setting for testing astroABC
	datacov - data covariance
	dfunc - distance function
	'''
	if  npart ==1 or npart < nparam:
		print ("\t Too few particles requested for the number of parameters. npart=%d nparam=%d") % (npart,nparam)
                print "\t exiting..."
                sys.exit(0)
	if nparam != len(priors):
		print ("\t Incorrect number of priors given for %d params") % nparam
		print "\t exiting..."
		sys.exit(0)
	if t1 < t2:
		print ("\t Tolerance levels should be [max,min] where max > min")
		print "\t exiting..."
		sys.exit(0)
	if dtype == "chi2" and datacov == None:
		print ("\t Data covariance matrix must be given for chi2 metric")
		print "\t exiting..."
		sys.exit(0)
	if dtype == "user" and not(dfunc):
		print ("\t A distance function must be specified, currently dfunc= %s") % dfunc
		print "\t exiting..."
		sys.exit(0)

def print_header(npart,niter,tol_type,tlevels,priors):
	'''print header to screen
	Input:
	npart - number of particles
	niter - number of iterations
	tol_type - decreasing tolerance level setting
	tlevels - max/min threshold levels
	priors - priors for parameters 
	'''
	print "\t \t"
	print "\t ########################     astroABC     ########################\t"
	print "\t \t"
	print "\t Npart=%d \t numt=%d \t tol=[%.4f,%.4f] %s"%(npart,niter,tlevels[0],tlevels[1], tol_type)
        print "\t Priors=", priors
	



