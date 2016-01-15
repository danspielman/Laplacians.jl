[TOC]

# Laplacians





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




