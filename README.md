# Laplacians.jl 

[![Build Status](https://travis-ci.org/danspielman/Laplacians.jl.svg?branch=master)](https://travis-ci.org/danspielman/Laplacians.jl)

Laplacians is a package containing graph algorithms, with an emphasis on tasks related to spectral and algebraic graph theory. It contains (and will contain more) code for solving systems of linear equations in graph Laplacians, low stretch spanning trees, sparsifiation, clustering, local clustering, and optimization on graphs.

All graphs are represented by sparse adjacency matrices. This is both for speed, and because our main concerns are algebraic tasks. It does not handle dynamic graphs. It would be very slow to implement dynamic graphs this way.

The documentation may be found in
[http://danspielman.github.io/Laplacians.jl/about/index.html](http://danspielman.github.io/Laplacians.jl/about/index.html).

This includes instructions for installing Julia, and some tips for how to start using it.  It also includes guidelines for Dan Spielman's collaborators.

For some examples of some of the things you can do with Laplacians, look at 

*  [this Julia notebook](http://github.com/danspielman/Laplacians.jl/blob/master/notebooks/FirstNotebook.ipynb).
*  [Low Stretch Spanning Trees](http://danspielman.github.io/Laplacians.jl/LSST/index.html
*  [Information about solving Laplacian equations](http://danspielman.github.io/Laplacians.jl/solvers/index.html)
*  And, try the chimera and wtedChimera graph generators.  They are designed to generate a wide variety of graphs so as to exercise code.
