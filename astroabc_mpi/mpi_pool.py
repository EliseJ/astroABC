try:
	from mpi4py import MPI
	MPI = MPI
except ImportError:
	MPI = None

class MpiPool(object):
	def __init__(self):
		self.comm = MPI.COMM_WORLD
		self.rank = MPI.COMM_WORLD.Get_rank()
		self.size = MPI.COMM_WORLD.Get_size()
		#self.mapFunction = mapFunction
		#print "\n Starting", self.rank,self.size

	def map(self, function, jobs, callback=None):
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
			 	#print "\n Worker_id", worker_id	
				req = self.comm.isend(job, dest=worker_id, tag=i)
				requests.append(req)

			MPI.Request.waitall(requests)

			results = []
			for i in range(njobs):
				worker_id = i % self.size + 1
				result = self.comm.recv(source=worker_id, tag=i)
				#print "\n Master recevied:", result, "from", worker_id


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
				#print "\n Master recevied:", result, "from", worker_id
				results[job] = result

                		if njobs_done < njobs:
					#print "\n NOT DONE< SENDING MORE JOBS..."
					job = jobs[njobs_done]
					#print "\n Sedning job", job
					i = njobs_done
					self.comm.isend(job, dest=worker_id, tag=i)
					njobs_done += 1


            		return results
			

	def worker(self):
		status = MPI.Status()
		while True:
			job = self.comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
			#print "\n \t\t recived",self.rank, job
			if job == -999:
				break
			if  isinstance(job, _func_wrapper):
				self.function = job.function
				continue

			result = self.function(job)
			#print "\n \t \t carried out work",self.rank, result
			self.comm.isend(result, dest=0, tag=status.tag)

	def close(self):
        	if self.rank == 0:
			requests=[]
            		for i in range(1,self.size):
                		req=self.comm.isend(-999, dest=i)
				requests.append(req)
                        MPI.Request.waitall(requests)

class _func_wrapper(object):
	def __init__(self, function):
		self.function = function
