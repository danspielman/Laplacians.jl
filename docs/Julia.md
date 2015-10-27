[TOC]

# Using Julia

To use the julia notebooks, you will need ipython and the IJulia package.  You should also get `jupyter`, which you should be able to install from ipython.
To install IJulia, you type `Pkg.add("IJulia")` from Julia.
Then, you just need to type `using IJulia` once.  This will tell jupyter about the Julia kernel.  To run Julia notebooks, you now type `jupyter notebook`.  You can select the kernel your new notebook is using.

Dan recommends installing the anaconda distribution of python.
You will then need to install some things from that, like
conda install jupyter (the new notebooks package)
conda install mathjax
conda install matplotlib 

This repository contains projects implemented in Julia by Dan Spielman's group.  While organizing by language is strange, we are trying it to help us learn the language.

There are two projects in here so far.  One, yinsGraph, is for doing graph theory and solving Laplacian systems.  It's documentation is at [yinsGraph/yinsGraph.md](yinsGraph/yinsGraph.md)

I would like to use the documentation pages in this root directory to discuss issues with how to get Julia to work well.  For now, I'll just ask some questions and write a little that I've figured out.  If you figure some out, please write it here.

## Converting to Julia 0.4

A description of the changes made in Julia 0.4 appears to be here
[https://github.com/JuliaLang/julia/blob/release-0.4/NEWS.md]
(https://github.com/JuliaLang/julia/blob/release-0.4/NEWS.md)

The first problem you will encounter when using Julia 0.4 is that it is stingier about its load path.  It won't load something from the current directory unless it is in the load path.  You can add a directory to the load path like

~~~julia
push!(LOAD_PATH,".")
~~~

or

~~~julia
push!(LOAD_PATH,"/Users/[your_username]/git/julia/yinsGraph")
~~~

To overcome this issue, you can add the above line to the end of the julia.rc file, found in 

~~~julia
/Applications/Julia-0.4.0-rc4.app/Contents/Resources/julia/etc/julia
~~~

Julia 0.4 lets you take advantage of docstrings.
For example, `?ringGraph` produces

~~~
The simple ring on n vertices
~~~

When having a multiline comment, make sure that lines don't have starting and trailing spaces. 
This will mess up the indentation when calling '?func_name'.

### Small details

* Julia 0.4 is trying to wean you off Matlab-like notation.  You should no longer create vectors like `[1:n]`.  Instead, you should type `collect(1:n)`

## Julia Notebooks
To get the Julia notebooks working, I presently type `jupyter notebook`.
I then select the kernel to be Julia-0.3.11.
It seems important to run this command from a directory that contains all the directories
that have notebooks that you will use.  In particular, I advise against "uploading" notebooks
from other directories.  That has only given me trouble.

The calico extensions that seem to be hosted at Brynmawr seem interesting.
I haven't yet figured out how to get them to work.
Here are the relevent links:

* http://jupyter.cs.brynmawr.edu/hub/dblank/public/Jupyter%20Help.ipynb
* http://jupyter.cs.brynmawr.edu/hub/dblank/public/Jupyter%20Notebook%20Users%20Manual.ipynb

To turn a notebook into html, you type something like

~~~
ipython nbconvert Laplacians.ipynb 
~~~

## Workflows

Julia has an IDE called Juno.  Both Dan and Serban have encountered some trouble with it: we have both found that it sometimes refuses to reload .jl code that we have written.  Please document workflows that you have found useful here:

#### Dan's current workflow:
* I use emacs (which has a mode for Julia) and the notebooks.
* I develop Julia code in a "temporary" file with a name like develX.jl.  While I am developing, this code is not included by the module to which it will eventually belong.
* After modifying code, I reload it with `include("develX.jl")`.  This works fine for reloading methods.  It is not a good way to reload modules or types.  So, I usually put the types either in a separate file, or in my julia notebook.

* I am writing this documention in MacDown.

#### Add your current workflow here:

## Things to be careful of (common bugs)

* Julia passes vectors and matrices to routines by reference, rather than by copying them.  If you type `x = y` when x and y are arrays, then this will make x a pointer to y.  If you want x to be a copy of y, type `x = copy(y)`.  This can really mess up matlab programmers.  I wrote many functions that were modifying their arguments without realizing it.

* On the other hand, if you type `x = x + y`, then x becomes a newly allocated vector and no longer refers to the original.  This is true even if you type `x += y`.  Here is an example that shows two of the possible behaviors, and the difference between what happens inside functions.

~~~julia

"Adds b in to a"
function addA2B(a,b)
    for i in 1:length(a)
        a[i] += b[i]
    end
end

"Fails to add b in to a"
function addA2Bfail(a,b)
	a += b
end

a = [1 0]
b = [2 2]
addA2B(a,b)
a

1x2 Array{Int64,2}:
 3  2
 
a = [1 0]
b = [2 2]
addA2Bfail(a,b)
a

1x2 Array{Int64,2}:
 1  0
 
a += b
a

1x2 Array{Int64,2}:
 3  2

~~~

* If you are used to programming in Matlab, you might be tempted to type a line like `for i in 1:10,`.  _Do not put extra commas in Julia!_  It will cause bad things to happen.

* Julia sparse matrix entries dissapear if they are set to 0. In order to overcome this, use the `setValue` function. `setValue(G, u, i, 0)` will set `weighti(G, u, i)` to 0 while also leaving `(u, nbri(G, u, i))` in the matrix.

## Useful Julia functions

I am going to make a short list of Julia functions/features that I find useful.  Please add those that you use often as well.

* docstrings: in the above example, I used a docstring to document each function.  You can get these by typing `?addA2B`.  You can also  [write longer docstrings and use markdown](http://julia.readthedocs.org/en/latest/manual/documentation/).  I suggest putting them in front of every function.  

* `methods(foo)` lists all methods with the name foo.
* `fieldnames(footype)` tells you all the fields of footype.  Note that this is 0.4.  In 0.3.11, you type `names(footype)`

~~~julia
julia> a = sparse(rand(3,3));
julia> fieldnames(a)
5-element Array{Symbol,1}:
 :m     
 :n     
 :colptr
 :rowval
 :nzval 
~~~


## Optimizing code in Julia

The best way that I've found of figuring out what's slowing down my code has been to use `@code_warntype`.  It only exists in version 4 of Julia.  For this reason, I keep one of those around.

Note that the first time you run a piece of code in Julia, it gets compiled.  So, you should run it on a small example before trying to time it.  Then, use `@time` to time your code.

I recommend reading the Performance Tips in the Julia documentation, not that I've understood all of it yet.

### Vectorization is Bad.
Julia is the anti-matlab in that vectorization is slow.
Still it is a good way to write your code the first time.
Here are some examples of code that adds one vector into another.
The first is vectorized, the second turns that into a loop, and the fastest uses BLAS.  Note that this was done in Julia 0.3.11.  The vectorized code is much faster, but still not fast, in 0.4.
_Also note that you have to run each routine once before it will be fast.  This is because it compiles it the first time your run it_

~~~julia
n = 10^7
a = rand(n)
b = rand(n)
@time a += b;

elapsed time: 0.155296017 seconds (80000312 bytes allocated)

a = rand(n)
b = rand(n)
@time add2(a,b);

elapsed time: 0.021190554 seconds (80 bytes allocated)

a = rand(n)
b = rand(n)
@time BLAS.axpy!(1.0,b,a);

elapsed time: 0.015894922 seconds (80 bytes allocated)

~~~

One reason that `a += b` was slow was that it seems to allocate a lot of memory.

## How should Julia packages be organized?

In yinsGraph, I decided to just make one big module called yinsGraph.jl.  It then includes a bunch of individual files, most of which contain many functions and types.  I think this is much nicer than making one file per functions, as some functions are very short.

I put the export statements in the main module.  The reason for this is that while developing code in a file, I don't include that in the module.  This way I can reload it as I change it without having to restart the kernel.  This does not seem to work as well for types.  I'm not sure why.

## How should notebooks play with Git?

The great thing about the notebooks is that they contain live code, so that you can play with them.  But, sometimes you get a version that serves as great documentation, and you don't want to klobber it my mistake later (or evern worse, have someone else klobber it).  Presumably if someone accidently commits a messed up version we can unwind that.  But, is there a good way to keep track of this?


