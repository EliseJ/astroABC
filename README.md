### **astroABC** ###

[![Build Status](https://travis-ci.com/EliseJ/astroABC.svg?token=LXdoQwTqixvxvKudVHQ7&branch=master)](https://travis-ci.com/EliseJ/astroABC)
[![Open Source Love](https://badges.frapsoft.com/os/mit/mit.svg?v=102)](https://github.com/EliseJ/astroABC/blob/master/LICENSE.txt)
 [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)



Version: 1.0.0

Author: Elise Jennings

### Description ###

astroABC is a Python implementation of an Approximate Bayesian Computation Sequential Monte Carlo (ABC SMC) sampler for parameter estimation. 

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

### Google group ###

Please post topics or suggestions to https://groups.google.com/forum/#!aboutgroup/astroabc 
astroabc@googlegroups.com 

### Wiki ###

For more information please read the [wiki](https://github.com/EliseJ/astroABC/wiki).

### Requirements ###

* numpy
* scipy
* mpi4py
* multiprocessing
* sklearn


### License ###

Copyright 2016 Elise Jennings

astroABC is free software made available under the MIT License. For details see the LICENSE.txt file.
