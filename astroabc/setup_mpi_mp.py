from mpi_pool import MpiPool
import multiprocessing as mp
import numpy as np
import copy_reg
import types

try:
	from mpi4py import MPI
	MPI = MPI
except ImportError:
	MPI = None


class Parallel(object):
	'''Class which sets up either an mpi or mp environment
	Input:
		mpi: Boolean True/False
		mp: Boolean True/False
		mpi_splitcomm: Boolean True/False, to split mpi comm
		num_proc: if mp, this is number of processors to use for mp.pool
		num_abc: if mpi_splitcomm, this is number of processors to use for abc sampling
		verbose: verbose level
	'''
	def __init__(self,mpi,mp,mpi_splitcomm,num_proc,num_abc,verbose):
		self.mpi=mpi
		self.mp=mp
		self.mpi_splitcomm=mpi_splitcomm
		self.num_proc=num_proc
		self.num_abc=num_abc
		self.verbose = verbose
		if self.mpi or self.mp or self.mpi_splitcomm:
			self.setup_mpi_mp()
		else:
			self.master=0


	def setup_mpi_mp(self):
		'''Method to setup pools and comms'''
		if self.mpi and not self.mpi_splitcomm:
			self.pool = MpiPool()
			self.master = self.pool.rank
			if self.pool.size ==1:
				print "\n \t Please run using mpirun -np #proc if you want to run using mpi. \n\texiting ..."
				sys.exit(0)
		elif self.mp:
			self.pool = mp.Pool(self.num_proc)
			self.master  = 0
		elif self.mpi_splitcomm:
			self.rank = MPI.COMM_WORLD.Get_rank()
			self.master = self.rank
			try:
				self.comm = MPI.COMM_WORLD
			except:
				print "\n \t Please run using mpirun -np #proc if you want to run using mpi. \n\texiting ..."
				sys.exit(0)

			self.group=self.comm.Get_group()

			if MPI.COMM_WORLD.Get_size() < self.num_abc:
				raise RuntimeError("Please specify at least -np X nodes for X particles. exiting... num_abc= %d" % num_abc)
				sys.exit(1)
			#set up group for abc sampler first
			self.abc_ranks=np.arange(0,self.num_abc)
			abc_group = self.group.Incl(self.abc_ranks)
			self.abc_comm = self.comm.Create(abc_group) 
			if self.rank in self.abc_ranks: #set up pool for abc sampler
				self.pool = MpiPool(comm=self.abc_comm)

			num_sim = MPI.COMM_WORLD.Get_size()-self.num_abc
			if num_sim ==0:
				raise RuntimeError("To the MPI comm please specify more nodes. Exiting... num_sim %d" % num_sim)
				sys.exit(1)
			if num_sim % (self.num_abc-1):   #every worker node doing abc has same number of node for simulation
				raise RuntimeError("Please specify equal num  mpi processors per node running abc in parallel. \
				exiting... num_sim % (num_abc-1) = %d" % (num_sim % (self.num_abc-1)))
				sys.exit(1)

			#set up groups for simulation pool
			self.sim_comms={} 
			self.all_sim_ranks = [[0]]
			for i in range(1,self.num_abc):
				start_rank=self.num_abc + int(num_sim/float(self.num_abc-1))*(i-1)
				end_rank=self.num_abc + int(num_sim/float(self.num_abc-1))*(i)
				sim_ranks=[i] 
				for j in np.arange(start_rank,end_rank):
					sim_ranks.append(j)
				sim_group = self.group.Incl(sim_ranks)
				self.sim_comms[i] = self.comm.Create(sim_group) 
				self.all_sim_ranks.append(sim_ranks)
			if self.verbose and self.master ==0:
				print "\n \t  MPI ranks for divided up:", self.all_sim_ranks
        

			self.sim_pool={}	
			for i,rank_list in enumerate(self.all_sim_ranks):
				if self.rank in rank_list and self.rank !=0:
					self.sim_pool = MpiPool(comm=self.sim_comms[i])
					
			if self.rank ==0: self.sim_pool=None #only node=0 won't have sim_pool


#http://stackoverflow.com/questions/25382455/python-notimplementederror-pool-objects-cannot-be-passed-between-processes
        def __getstate__(self):
                self_dict = self.__dict__.copy()
                del self_dict['pool']
                return self_dict

        def __setstate__(self, state):
                self.__dict__.update(state)


def _pickle_method(method):
        """
        http://stackoverflow.com/questions/11726809/
        python-efficient-workaround-for-multiprocessing-a-function-that-is-a-data-member
        """
        func_name = method.im_func.__name__
        obj = method.im_self
        cls = method.im_class
        cls_name = ''
        if func_name.startswith('__') and not func_name.endswith('__'):
                cls_name = cls.__name__.lstrip('_')
        if cls_name:
                func_name = '_' + cls_name + func_name
        return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
        for cls in cls.mro():
                try:
                        func = cls.__dict__[func_name]
                except KeyError:
                        pass
                else:
                        break
        return func.__get__(obj, cls)

copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


