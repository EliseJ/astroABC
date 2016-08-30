
<a href="https://github.com/EliseJ/astroABC"><img src="https://github.com/EliseJ/astroABC/blob/master/abc_logo.001.jpeg"
align="left" hspace="10" vspace="6"></a>

**astroABC**

<br>
[![Build Status](https://travis-ci.org/EliseJ/astroABC.svg?branch=master)](https://travis-ci.com/EliseJ/astroABC)
[![Latest Version](http://img.shields.io/pypi/v/astroabc.svg?style=flat)](https://pypi.python.org/pypi/astroabc/)
[![Open Source Love](https://badges.frapsoft.com/os/mit/mit.svg?v=102)](https://github.com/EliseJ/astroABC/blob/master/LICENSE.txt)
 [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/EliseJ/astroABC/issues)



Author: Elise Jennings

[arxiv:1608.07606](http://arxiv.org/abs/1608.07606)

**astroABC** is a Python implementation of an Approximate Bayesian Computation Sequential Monte Carlo (ABC SMC) sampler for parameter estimation. 


<br>
<br>

## Key features ##

- Parallel sampling using MPI or multiprocessing
- MPI communicator can be split so both the sampler, and simulation launched by each particle, can run in parallel
- A Sequential Monte Carlo sampler (see e.g. [Toni et al. 2009], [Beaumont et al. 2009], [Sisson & Fan 2010])
[Toni et al. 2009]:https://arxiv.org/abs/0901.1925
[Sisson & Fan 2010]:http://arxiv.org/abs/1001.2058
[Beaumont et al. 2009]:https://arxiv.org/abs/0805.2256
- A method for iterative adapting tolerance levels using the qth quantile of the distance for t iterations ([Turner & Van Zandt (2012)])
- Scikit-learn covariance matrix estimation using [Ledoit-Wolf shrinkage] for singular matrices
[Ledoit-Wolf shrinkage]:http://scikit-learn.org/stable/modules/covariance.html
- A module for specifying particle covariance using method proposed by [Turner & Van Zandt (2012)], optimal covariance matrix  for a multivariate normal perturbation kernel, local covariance estimate using scikit-learn  KDTree method for nearest neighbours ([Filippi et al 2013]) and a weighted covariance (Beaumont et al 2009)
[Turner & Van Zandt (2012)]:http://link.springer.com/article/10.1007/s11336-013-9381-x
[Filippi et al 2013]:https://arxiv.org/abs/1106.6280
- Restart files output frequently so an interrupted run can be resumed at any iteration
- Output and restart files are backed up every iteration 
- User defined distance metric and simulation methods
- A class for specifying heterogeneous parameter priors 
- Methods for drawing from any non-standard prior PDF e.g using Planck/WMAP chains 
- A module for specifying a constant, linear, log or exponential tolerance level
- Well-documented examples and sample scripts


### Wiki ###

For more information please read the [wiki](https://github.com/EliseJ/astroABC/wiki).

### Installing ###

Install astroABC using pip

```
$ pip install astroabc
```

or git clone the repository using the url above. 
Check the dependencies listed in the next section are installed.

### Dependencies ###

* numpy
* scipy
* mpi4py
* multiprocessing
* sklearn

Python distributions like [Anaconda] have most of what is needed. 
You can then conda install or pip install all of the required dependencies.

```
$ conda install  numpy scipy scikit-learn mpi4py
$ pip install numpy scipy scikit-learn mpi4py
```

[Anaconda]:https://www.continuum.io/downloads

### License ###

Copyright 2016 Elise Jennings

astroABC is free software made available under the MIT License. For details see the LICENSE.txt file.
