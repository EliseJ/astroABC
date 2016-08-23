try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser  # ver. < 3.0
import numpy as np

config = ConfigParser()

class ParseParams():
    def __init__(self,inifile):
        self.inifile =inifile
        config.read(inifile)
        self.get_datafiles()
        self.get_abcparams()
        self.get_priors()

    def get_datafiles(self):
        #self.dirname = config.get('datafiles', 'dirname')
        self.datafile = config.get('datafiles', 'data_file')
        self.model_type = config.get('datafiles', 'model_type')
        self.true_p0 = np.float(config.get('datafiles', 'true_param_0'))
        self.true_p1 = np.float(config.get('datafiles', 'true_param_1'))
        self.nsamples = config.getint('datafiles', 'nsamples')

    def get_abcparams(self):
        self.nparam = config.getint('abc', 'nparam')
        self.npart = config.getint('abc', 'npart')
        self.niter = config.getint('abc', 'niter')
        self.tlevels = [np.float(k) for k in config.get('abc', 'tlevels').split(',')]
        self.tol_type = config.get('abc', 'tol_type')
        self.verbose = config.getint('abc', 'verbose')
        self.adapt_t = config.getboolean('abc', 'adapt_t')
        self.threshold = config.getint('abc', 'threshold')
        self.pert_kernel= config.getint('abc', 'pert_kernel')
        self.var_method= config.getint('abc', 'variance_method')
        self.k_near= config.getint('abc', 'k_near')
        self.dist_type= config.get('abc', 'dist_type')
        self.mpi = config.getboolean('abc', 'mpi')
        self.mpi_splitcomm = config.getboolean('abc', 'mpi_splitcomm')
        self.num_abc= config.getint('abc', 'num_abc')
        self.mp = config.getboolean('abc', 'mp')
        self.num_proc= config.getint('abc', 'num_proc')
        self.outfile= config.get('abc', 'outfile')
        self.restart= config.get('abc', 'restart')
        self.from_restart= config.getboolean('abc', 'from_restart')

    def get_priors(self):
        p = []
        hyp = []
        pnames=[]
        for i in range(self.nparam):
            p.append(config.get('priors', 'priorname'+str(i)))
            temp = [np.float(j) for j in config.get('priors', 'hyperp'+str(i)).split(",")]
            hyp.append(temp)
            pnames.append(config.get('priors', 'param_name'+str(i)))
        self.priors = p
        self.hyperp = hyp
        self.pnames = pnames
            

