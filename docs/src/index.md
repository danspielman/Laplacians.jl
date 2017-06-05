# Laplacians.jl 

[![Build Status](https://travis-ci.org/danspielman/Laplacians.jl.svg?branch=master)](https://travis-ci.org/danspielman/Laplacians.jl)
[![codecov](https://codecov.io/gh/danspielman/Laplacians.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/danspielman/Laplacians.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://danspielman.github.io/Laplacians.jl/latest)



Laplacians is a package containing graph algorithms, with an emphasis on tasks related to spectral and algebraic graph theory. It contains (and will contain more) code for solving systems of linear equations in graph Laplacians, low stretch spanning trees, sparsifiation, clustering, local clustering, and optimization on graphs.

All graphs are represented by sparse adjacency matrices. This is both for speed, and because our main concerns are algebraic tasks. It does not handle dynamic graphs. It would be very slow to implement dynamic graphs this way.

The documentation may be found by clicking on one of the "docs" links above.


This includes instructions for installing Julia, and some tips for how to start using it.  It also includes guidelines for Dan Spielman's collaborators.

For some examples of some of the things you can do with Laplacians, look at 

*  [this Julia notebook](http://github.com/danspielman/Laplacians.jl/blob/master/notebooks/FirstNotebook.ipynb).
*  [Low Stretch Spanning Trees](LSST.md)
*  [Information about solving Laplacian equations](usingSolvers.md)
*  [An example of sparsification](http://github.com/danspielman/Laplacians.jl/blob/master/notebooks/Sparsification Demo.ipynb)
*  And, try the chimera and wtedChimera graph generators.  They are designed to generate a wide variety of graphs so as to exercise code.

If you want to solve Laplacian equations, we recommend the KMPLapSolver.  For SDD equations, we recommend the KMPSDDSolver.

The algorithms provide by Laplacians.jl include:

* `akpw`, a heuristic for computing low stretch spanning trees written by Daniel Spielman, inspired by the algorithm from the paper "A graph-theoretic game and its application to the k-server problem" by Alon, Karp, Peleg, and West, <i>SIAM Journal on Computing</i>, 1995.
* `approxCholLap`: a fast heuristic for solving Laplacians equations written by Daniel Spielman, based on the paper "Approximate Gaussian Elimination for Laplacians: Fast, Sparse, and Simple" by Rasmus Kyng and Sushant Sachdeva, FOCS 2016. 
* `sparsify`, an implementation of sparsification by effective resistance sampling, following Spielman and Srivastava.
* `KMPLapSolver` and `KMPSDDSolver`: linear equation solvers based on the paper "Approaching optimality for solving SDD systems" by Koutis, Miller, and Peng, <i>SIAM Journal on Computing</i>, 2014.
* `samplingSDDSolver` and `samplingLapSolver`, based on the paper "Approximate Gaussian Elimination for Laplacians: Fast, Sparse, and Simple" by Rasmus Kyng and Sushant Sachdeva, FOCS 2016. 
* `chimera` and `wtedChimera` graph generators for testing graph algorithms, by Daniel Spielman.
* Local Graph Clustering Heuristics, implemented by Serban Stan, including `prn` a version of PageRank Nibble based on "Using PageRank to Locally Partition a Graph", <i>Internet Mathematics</i> and `LocalImprove` based on "Flow-Based Algorithms for Local Graph Clustering" by Zeyuan Allen-Zhu and Lorenzo Orecchia, SODA 2014.

