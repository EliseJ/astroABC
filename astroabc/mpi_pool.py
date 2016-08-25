try:
	from mpi4py import MPI
	MPI = MPI
except ImportError:
	MPI = None

class MpiPool(object):
	'''mpi pool class'''
	def __init__(self,comm=None):
		'''initialize a communicator and set the rank and size '''
		if comm==None:
			self.comm = MPI.COMM_WORLD
		else:
			self.comm=comm
		self.rank = self.comm.Get_rank()
		self.size = self.comm.Get_size()
		self.status = MPI.Status()

	def map(self, function, jobs):
		'''
		Map function to perform similar function to multiprocessing mp.pool map()
		Input:
		function  - function to be mapped
		jobs - array of jobs to be assigned to each proc
		
		'''
		njobs = len(jobs)
		self.function = function
		
		
		# If not the master just wait for instructions.
		if not self.rank == 0:
			self.worker()
			return 

		F = _func_wrapper(function)
		req = [self.comm.isend(F, dest=i) for i in range(1,self.size)]
		MPI.Request.waitall(req)
		
		if  (njobs <= self.size-1):
			requests = []
			for i, job in enumerate(jobs):
               			worker_id = i % self.size + 1
				req = self.comm.isend(job, dest=worker_id, tag=i)
				requests.append(req)

			MPI.Request.waitall(requests)

			results = []
			for i in range(njobs):
				worker_id = i % self.size + 1
				result = self.comm.recv(source=worker_id, tag=i)
				results.append(result)
			return results

		else:
			for i in range(1, self.size):
				req = self.comm.isend(jobs[i-1], dest=i, tag=i)

			njobs_done = self.size-1 # don't use master for function
			results = [None]*njobs
			for job in range(njobs):
				status = MPI.Status()
				result = self.comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG, status=status)
				worker_id = status.source
				results[job] = result

                		if njobs_done < njobs:
					job = jobs[njobs_done]
					i = njobs_done
					self.comm.isend(job, dest=worker_id, tag=i)
					njobs_done += 1


            		return results
			

	def worker(self):
		'''worker processor method which waits for map data or end signal from master'''
		status = MPI.Status()
		while True:
			job = self.comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
			if job == -999: #end signal
				break
			if  isinstance(job, _func_wrapper):
				self.function = job.function
				continue

			result = self.function(job)
			self.comm.isend(result, dest=0, tag=status.tag)

	def close(self):
		'''method to close the pool, master sends end signal'''
        	if self.rank == 0:
			requests=[]
            		for i in range(1,self.size):
                		req=self.comm.isend(-999, dest=i)
				requests.append(req)
                        MPI.Request.waitall(requests)

class _func_wrapper(object):
	'''wrapper class for function
	modified from cosmosis / cosmosis / runtime / mpi_pool.py
	'''
	def __init__(self, function):
		self.function = function
