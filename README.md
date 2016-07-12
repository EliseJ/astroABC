### **astroABC** ###
## author: Elise Jennings ##
## elise@fnal.gov ##


### Description ###

astroABC is a Python implementation of an Approximate Bayesian Computation Sequential Monte Carlo (ABC SMC) sampler for parameter estimation. 

astroABC is applicable to a large suite of problems. 

### Key features include ###

- Parallel sampling using MPI or multiprocessing
- A Sequential Monte Carlo sampler following Beaumont et al. 2009
- A method for iterative adapting tolerance levels using the qth quantile of the distance for t iterations (Drovandi & Pettitt 2011)
- Scikit-learn covariance matrix estimation using Ledoit-Wolf shrinkage for singular matrices
- A module for specifying particle covariance using method proposed by Turner & Van Zandt (2012), optimal covariance matrix  for a multivariate normal perturbation kernel (Filippi et al 2012) and a weighted covariance metric (Beaumont et al 2009)
- Restart files output frequently so an interrupted sampling run can be resumed at any iteration
- Output and restart files are backed up every iteration before new output is written
- User defined distance metric and simulation methods
- A module for specifying heterogeneous parameter priors 
- A module for specifying a constant, linear, log or exponential tolerance level
- Well-documented examples and sample scripts



For more information please read the [wiki](https://bitbucket.org/elisejennings/astroabc_mpi/wiki/Home).

### Requirements ###

* mpi4py
* multiprocessing
* sklearn
* scipy
* numpy


### License ###

Copyright 2016 Elise Jennings

astroABC is free software made available under the MIT License. For details see the LICENSE.txt file.
