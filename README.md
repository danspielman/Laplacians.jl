# Laplacians.jl 

[![Build Status](https://travis-ci.org/danspielman/Laplacians.jl.svg?branch=master)](https://travis-ci.org/danspielman/Laplacians.jl)
[![codecov](https://codecov.io/gh/danspielman/Laplacians.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/danspielman/Laplacians.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://danspielman.github.io/Laplacians.jl/latest)



Laplacians is a package containing graph algorithms, with an emphasis on tasks related to spectral and algebraic graph theory. It contains (and will contain more) code for solving systems of linear equations in graph Laplacians, low stretch spanning trees, sparsifiation, clustering, local clustering, and optimization on graphs.

All graphs are represented by sparse adjacency matrices. This is both for speed, and because our main concerns are algebraic tasks. It does not handle dynamic graphs. It would be very slow to implement dynamic graphs this way.

The documentation may be found by clicking on one of the "docs" links above.


## Current Development Version

To get the current version of the master branch, run `Pkg.checkout("Laplacians")`

## Version 0.2.2, December 28, 2017

Fixed two bugs: one in shortestPaths, and one that prevented passing some parameters to approxCholSddm.  Improved the documentation for solving linear equations.

## Version 0.2.1, September 18, 2017

Fixed a bug in `approxCholSddm` that caused it to be slow.

## Version 0.2.0, July 17, 2017

This version is compatible with Julia 0.6.  It will not work with
Julia 0.5.X.

Changes:

* Added `approxCholSddm`, a wrapper of `approxCholLap` that solves
  SDDM systems.

## Version 0.1.4, June 6, 2017

This is the current version.  It is what you retrieve when you run `Pkg.add("Laplacians")`. 

Changes:

* Added `sparsify`, an implementation of sparsification by effective resistance sampling, following Spielman and Srivastava.
* Added `approxQual` and `conditionNumber` for checking how well one graph approximates another.
* Fixed a bug in the solution of Laplacian systems in disconnected graphs.

## Version 0.1.3, June 2, 2017

Major Changes:

* Changed the name of the approximate Cholesky solver from `edgeElimLap` to `approxCholLap`.  Made improvements in this solver.
* Improved PCG so that it can now detect stagnation.  Made options to do this even better when using it with a good preconditioner, like `approxCholLap`.
* Added in code for comparing the running times of solvers.  The difficulty here is that we need to stop them if they run too long.  Added code to do this with threads inside Julia, and with `gtimeout` when calling Matlab to use icc, CMG, or LAMG.

## Version 0.1.2, April 2, 2017


Major Changes:

* added `edgeElimLap` - a fast Laplacian solver.
* fixed a bug in the unweighted version of `akpw`. 

## Version 0.1.1, December 26, 2016

Changelist:

* All of the linear equation solvers now have the same interface, and the Laplacian solvers work for disconnected graphs. 
* Some support for calling solvers from Matlab has been added.
* Documentation is now through Documenter.jl.

## Version 0.0.3 / 0.1.0, November 20, 2016

Versions 0.0.3 and 0.1.0 are the same.
These versions works with Julia 0.5.

Warning: the behavior of chimera and wtedChimera differs between Julia 0.4 and Julia 0.5 because randperm acts differently in these.

## Version 0.0.2, November 19, 2016

This is the version that works with Julia 0.4.
It was captured right before the upgrade to Julia 0.5



# Acknowledgements

The development of this package has been supported in part by the
National Science Foundation Award CCF-1562041  and by the Office of Naval
Research Award N00014-16-1-2374.
