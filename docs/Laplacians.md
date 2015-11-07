[TOC]

# Laplacians




## Installation

First, you will need Julia.
You will also need a number of Julia packages.  
You install these like

~~~julia
Pkg.add("PyPlot")
Pkg.add("DataStructures")
~~~

Once these are installed, you can use Laplacians.
Right now, Laplacians is a module, not a package.
So, we will need to do a little more to get it started.
In the directory where Laplacians resides, type the following:

~~~julia
push!(LOAD_PATH,"src")
using Laplacians
~~~

Instead of adding to the load path every time you use Julia, you could put the following line in a file called .juliarc.jl that lives in your home directory:

~~~julia
push!(LOAD_PATH,"[path_to_laplacians]/Laplacians.jl/src")
~~~

In my case, the path to laplacians is `/Users/spielman/git/`.
Then, when you want to use the module, you just need to type `using Laplacians`.

You may need to install matplotlib in python before PyPlot.
Look at this page for more information: https://github.com/stevengj/PyPlot.jl

If you discover that you need any other packages, please list them above.

Other recommended (but not necessary) packages are:

* Optim

To see if it is working, try something like this:

~~~julia
a = chimera(100,6)
spectralDrawing(a)
~~~

or

~~~julia
a = generalizedNecklace(grid2(6),grid2(3),2)
spectralDrawing(a)
~~~

## To use Laplacians

Examples of how to do many things in yinsGraph may be found in the IJulia notebooks.  These have the extensions .ipynb.  When they look nice, I think it makes sense to convert them to .html.

Right now, the notebooks worth looking at are:

* [yinsGraph](yinsGraph.html) - usage, demo, and speed tests (Laplacians was previously called yinsGraph)
* [Solvers](solvers.md) - code for solving equations.  How to use direct methods, conjugate gradient, and a preconditioned augmented spanning tree solver.


(I suggest that you open the html in your browser)




