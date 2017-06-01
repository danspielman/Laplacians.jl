# Laplacians.jl 

[![Build Status](https://travis-ci.org/danspielman/Laplacians.jl.svg?branch=master)](https://travis-ci.org/danspielman/Laplacians.jl)
[![codecov](https://codecov.io/gh/danspielman/Laplacians.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/danspielman/Laplacians.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://danspielman.github.io/Laplacians.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://danspielman.github.io/Laplacians.jl/latest)



Laplacians is a package containing graph algorithms, with an emphasis on tasks related to spectral and algebraic graph theory. It contains (and will contain more) code for solving systems of linear equations in graph Laplacians, low stretch spanning trees, sparsifiation, clustering, local clustering, and optimization on graphs.

All graphs are represented by sparse adjacency matrices. This is both for speed, and because our main concerns are algebraic tasks. It does not handle dynamic graphs. It would be very slow to implement dynamic graphs this way.

The documentation may be found by clicking on one of the "docs" links above.

