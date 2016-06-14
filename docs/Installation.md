[TOC]

# Installation

Before you can use Laplacians, you need Julia.
So, we'll begin with instructions for installing Julia.  I (Dan S.) found that it worked best if I installed Python first.  So, I'll suggest that you do the same.

All of these instruction assume you are using a Mac.  

## Python

Install python.  I recommend the anaconda distribution [https://www.continuum.io/downloads](https://www.continuum.io/downloads).

Once you install python, you are going to want two packages: a plotting package that Julia will use, and jupyter notebooks for interacting with Julia.  Install them like

~~~
conda install matplotlib
conda install mathjax
conda install jupyter
~~~

## Julia

You can get Julia from 
[http://julialang.org/](http://julialang.org/).  I haven't had much luck with Juno, so I recommend the plain Julia installation.

Once you have this, you will want Julia notebooks.  To install this, run `julia` and type

~~~julia
julia> Pkg.add("IJulia")
julia> using IJulia
~~~

This will install the package, and put the current julia kernel into [jupyter](http://jupyter.org/).  In the future, you can launch the Julia notebook by typing (in a terminal)

~~~ 
jupyter notebook
~~~

## Laplacians

In theory, all you need to do now is type either

~~~julia
julia> Pkg.clone("git://github.com/danspielman/Laplacians.jl.git")
~~~

Or, in a few days,

~~~julia
julia> Pkg.add("Laplacians")
~~~

This should add all the packages upon which Laplacians explicitly depends.  


## Using Laplacians

To actually use the Laplacians package in a Julia session, you must type

~~~julia
julia> using Laplacians 
~~~

To see if Laplacians is working, try typing

~~~julia
a = chimera(100,6)
spectralDrawing(a)
~~~

or

~~~julia
a = generalizedNecklace(grid2(6),grid2(3),2)
spectralDrawing(a)
~~~
