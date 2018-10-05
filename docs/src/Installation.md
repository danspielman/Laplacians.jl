# Installation

**Caveat: These instructions are old (circa Julia 0.5)  The same ideas apply, but they need updating.  You can probably find a better reference now.** 


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
[http://julialang.org/](http://julialang.org/).  
If you are using a Mac, you may wish to create a symnolic link to the Julia executable so that you can call it from a terminal.  For example, you can do this like:

~~~sh
cd /usr/local/bin/
ln -s julia /Applications/Julia-0.5.app/Contents/Resources/julia/bin/julia
~~~

Once you have this, you will want Julia notebooks.  To install this, run `julia` and type

~~~julia
julia> Pkg.add("IJulia")
julia> using IJulia
~~~

This will install the package, and put the current julia kernel into [jupyter](http://jupyter.org/).  In the future, you can launch the Julia notebook by typing (in a terminal)

~~~ 
jupyter notebook
~~~

Most users will also want to install PyPlot, if you did not already.
To do that, type

~~~julia
Pkg.add("PyPlot")
~~~

## Laplacians

In theory, all you need to do now is type either

~~~julia
julia> Pkg.add("Laplacians")
~~~

To use the package, you then type

~~~julia
julia> using Laplacians
~~~

The one catch is with the functions for drawing graphs.  These require PyPlot.  If you did not install it before typing `Pkg.add("PyPlot")`, then you can either install it now or disable the plotting routines in Laplacians.

If you do not want to load PyPlot, then either set the environment variable `LAPLACIANS_NOPLOT` to `true` in bash, like

~~~
$ export LAPLACIANS_NOPLOT=true
~~~

or, set the variable inside Julia, like 

~~~julia
julia> LAPLACIANS_NOPLOT = true
~~~

before typing `using Laplacians`.
Similarly, you can avoid loading PyAmg by setting

~~~julia
julia> LAPLACIANS_NOAMG = true
~~~

(note that these are not the same variable: the environment variable in Julia is available as `ENV["LAPLACIANS_NOPLOT"]`.
Actually, defining these variables to anything will have the same
effect.  So, setting them to false has the same effect as setting them
to true.



To see if Laplacians is working, try typing

~~~julia
a = chimera(100,6)
spectral_drawing(a)
~~~

or

~~~julia
a = generalizedNecklace(grid2(6),grid2(3),2)
spectral_drawing(a)
~~~
