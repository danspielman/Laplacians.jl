var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "About",
    "title": "About",
    "category": "page",
    "text": ""
},

{
    "location": "#Laplacians.jl-1",
    "page": "About",
    "title": "Laplacians.jl",
    "category": "section",
    "text": "(Image: Build Status) (Image: codecov) (Image: )Laplacians is a package containing graph algorithms, with an emphasis on tasks related to spectral and algebraic graph theory. It contains (and will contain more) code for solving systems of linear equations in graph Laplacians, low stretch spanning trees, sparsifiation, clustering, local clustering, and optimization on graphs.All graphs are represented by sparse adjacency matrices. This is both for speed, and because our main concerns are algebraic tasks. It does not handle dynamic graphs. It would be very slow to implement dynamic graphs this way.The documentation may be found by clicking the \"docs\" link above.This includes instructions for installing Julia, and some tips for how to start using it.  It also includes guidelines for Dan Spielman\'s collaborators.For some examples of some of the things you can do with Laplacians, look atthis Julia notebook.\nLow Stretch Spanning Trees\nInformation about solving Laplacian equations\nAn example of sparsification\nAnd, try the chimera and wtedChimera graph generators.  They are designed to generate a wide variety of graphs so as to exercise code.If you want to solve Laplacian equations, we recommend approxchol_lap.The algorithms provide by Laplacians.jl include:akpw, a heuristic for computing low stretch spanning trees written by Daniel Spielman, inspired by the algorithm from the paper \"A graph-theoretic game and its application to the k-server problem\" by Alon, Karp, Peleg, and West, <i>SIAM Journal on Computing</i>, 1995.\napproxchol_lap: a fast heuristic for solving Laplacians equations written by Daniel Spielman, based on the paper \"Approximate Gaussian Elimination for Laplacians: Fast, Sparse, and Simple\" by Rasmus Kyng and Sushant Sachdeva, FOCS 2016.   For SDDM systems, use approxchol_sddm.\nharmonic_interp: Harmonic Interpolation on graphs.  Minimizes the Laplacians quadratic form subject to fixing values at certain vertices.\nsparsify, an implementation of sparsification by effective resistance sampling, following Spielman and Srivastava.\nKMPLapSolver and KMPSDDSolver: linear equation solvers based on the paper \"Approaching optimality for solving SDD systems\" by Koutis, Miller, and Peng, <i>SIAM Journal on Computing</i>, 2014.\nsamplingSDDSolver and samplingLapSolver, based on the paper \"Approximate Gaussian Elimination for Laplacians: Fast, Sparse, and Simple\" by Rasmus Kyng and Sushant Sachdeva, FOCS 2016.\nchimera and wted_chimera graph generators for testing graph algorithms, by Daniel Spielman.\nLocal Graph Clustering Heuristics, implemented by Serban Stan, including prn a version of PageRank Nibble based on \"Using PageRank to Locally Partition a Graph\", <i>Internet Mathematics</i> and LocalImprove based on \"Flow-Based Algorithms for Local Graph Clustering\" by Zeyuan Allen-Zhu and Lorenzo Orecchia, SODA 2014."
},

{
    "location": "#Current-Development-Version-1",
    "page": "About",
    "title": "Current Development Version",
    "category": "section",
    "text": "To get the current version of the master branch, run pkg> add Laplacians#master"
},

{
    "location": "#Version-1.2.0-1",
    "page": "About",
    "title": "Version 1.2.0",
    "category": "section",
    "text": "This version is compatible with Julia 1.4 and 1.5, but not earlier versions.Changes:Added two graph generators: complete_bipartite_graph, star_graph.\nAdded a function line_graph that computes the line graph of an input graph."
},

{
    "location": "#Version-1.1.1-1",
    "page": "About",
    "title": "Version 1.1.1",
    "category": "section",
    "text": "Change: minor bug fix for spectral graph drawing.Verified compatibility with Julia 1.2."
},

{
    "location": "#Version-1.1.0-1",
    "page": "About",
    "title": "Version 1.1.0",
    "category": "section",
    "text": "Changes:Updating to use Julia\'s new Registrator.\nAdded harmonic_interp to interpolate harmonic functions on graphs.  This is the fundamental routine used in Label Propagation / Semi-Supervised Learning on Graphs.\nAdded a function read_graph to replace readIJ and readIJV.  It is a little more robust.\nCleaned up maxflow so that it now returns a flow and cut as a matrix and set.\nMade pcg a little more general.\nAdded fiedler for computing Fiedler vectors and values.  That is, the smallest nonzero eigenvalue of the Laplacian.\nFixed a bug in thicken that could cause it to loop forever, and cause chimera to do the same.\nChanged the graph drawing code to use Plots instead of PyPlot."
},

{
    "location": "#Version-1.0.1-1",
    "page": "About",
    "title": "Version 1.0.1",
    "category": "section",
    "text": "Changes:Added latin_square_graph and latin_square.\nAllow plot_graph to plot in 3D.\nFixed performance bug due to lazy matrix transpose.\nChanged more function names to agree with Julia naming conventions.\nUpdate documentation and examples."
},

{
    "location": "#Version-1.0.0-1",
    "page": "About",
    "title": "Version 1.0.0",
    "category": "section",
    "text": "This version works with Julia version 1.0.0."
},

{
    "location": "#Verson-0.3.1-1",
    "page": "About",
    "title": "Verson 0.3.1",
    "category": "section",
    "text": "Changes:The major change in this version is to the chimera and wted_chimera graph generators.  They are now faster, and incorporate two-lifts and thickening.  The old versions, using the pseudorandom generator from Julia V0.6 and Versions 0.2 of Laplacians, may be accessed by using the flag ver=Laplacians.V06, as in\na = chimera(2000, 1, ver=Laplacians.V06)\nThere do seem to be differences in the very low order bits of graphs generated by wted_chimera with the V06 option and those generated in Julia V0.6.  Not sure why.\nThe old generator is obtained by using the RandomV06 package for Julia.\nChanged the names of many functions to bring closer to the Julia standard naming scheme.  New names are emptygraph, pathgraph, ringgraph, completegraph, generalizedring, randgenring, productgraph, joingraphs, twolift ...  Set deprecation warnings for the old names.\nMoved lex.jl to the directory buggy, as on further testing we found bugs in it.\ndropped wGrid3, as it produced a 4d grid so probably wasn\'t being used anyway.  Dropped wGrid2 also."
},

{
    "location": "#Version-0.3.0,-July-18-(or-so),-2017-1",
    "page": "About",
    "title": "Version 0.3.0, July 18 (or so), 2017",
    "category": "section",
    "text": "This is the first version that is compatible with Julia 0.7.  Other changes:Dropped support for samplingSDDM and samplingLap solvers.\nThe behavior of rand in Julia 0.7 is different, and this has changed the behavior of chimera.  So, the chimera graphs generated in Version 0.3.0 and beyond will be different from those before."
},

{
    "location": "#Version-0.2.2,-December-28,-2017-1",
    "page": "About",
    "title": "Version 0.2.2, December 28, 2017",
    "category": "section",
    "text": "Fixed two bugs: one in shortestPaths, and one that prevented passing some parameters to approxchol_sddm.  Improved the documentation for solving linear equations."
},

{
    "location": "#Version-0.2.1,-September-18,-2017-1",
    "page": "About",
    "title": "Version 0.2.1, September 18, 2017",
    "category": "section",
    "text": "Fixed a bug in approxchol_sddm that caused it to be slow."
},

{
    "location": "#Version-0.2.0,-July-17,-2017-1",
    "page": "About",
    "title": "Version 0.2.0, July 17, 2017",
    "category": "section",
    "text": "This version is compatabile with Julia 0.6.  It will not work with Julia 0.5.X.Changes:Added approxchol_sddm, a wrapper of approxchol_lap that solves SDDM systems."
},

{
    "location": "#Version-0.1.4,-June-6,-2017-1",
    "page": "About",
    "title": "Version 0.1.4, June 6, 2017",
    "category": "section",
    "text": "This is the current version.  It is what you retrieve when you run Pkg.add(\"Laplacians\").Changes:Added sparsify, an implementation of sparsification by effective resistance sampling, following Spielman and Srivastava.\nAdded approxQual and conditionNumber for checking how well one graph approximates another.\nFixed a bug in the solution of Laplacian systems in disconnected graphs."
},

{
    "location": "#Version-0.1.3,-June-2,-2017-1",
    "page": "About",
    "title": "Version 0.1.3, June 2, 2017",
    "category": "section",
    "text": "Major Changes:Changed the name of the approximate Cholesky solver from edgeElimLap to approxchol_lap.  Made improvements in this solver.\nImproved PCG so that it can now detect stagnation.  Made options to do this even better when using it with a good preconditioner, like approxchol_lap.\nAdded in code for comparing the running times of solvers.  The difficulty here is that we need to stop them if they run too long.  Added code to do this with threads inside Julia, and with gtimeout when calling Matlab to use icc, CMG, or LAMG."
},

{
    "location": "#Version-0.1.2,-April-2,-2017-1",
    "page": "About",
    "title": "Version 0.1.2, April 2, 2017",
    "category": "section",
    "text": "This is the current version.  It is what you retrieve when you run Pkg.add(\"Laplacians\").Major Changes:added edgeElimLap - a fast Laplacian solver.\nfixed a bug in the unweighted version of akpw."
},

{
    "location": "#Version-0.1.1,-December-26,-2016-1",
    "page": "About",
    "title": "Version 0.1.1, December 26, 2016",
    "category": "section",
    "text": "Changelist:All of the linear equation solvers now have the same interface, and the Laplacian solvers work for disconnected graphs.\nSome support for calling solvers from Matlab has been added.\nDocumentation is now through Documenter.jl."
},

{
    "location": "#Version-0.0.3-/-0.1.0,-November-20,-2016-1",
    "page": "About",
    "title": "Version 0.0.3 / 0.1.0, November 20, 2016",
    "category": "section",
    "text": "Versions 0.0.3 and 0.1.0 are the same. These versions works with Julia 0.5.Warning: the behavior of chimera and wtedChimera differs between Julia 0.4 and Julia 0.5 because randperm acts differently in these."
},

{
    "location": "#Version-0.0.2,-November-19,-2016-1",
    "page": "About",
    "title": "Version 0.0.2, November 19, 2016",
    "category": "section",
    "text": "This is the version that works with Julia 0.4. It was captured right before the upgrade to Julia 0.5"
},

{
    "location": "Installation/#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "Installation/#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "Caveat: These instructions are old (circa Julia 0.5)  The same ideas apply, but they need updating.  You can probably find a better reference now. Before you can use Laplacians, you need Julia. So, we\'ll begin with instructions for installing Julia.  I (Dan S.) found that it worked best if I installed Python first.  So, I\'ll suggest that you do the same.All of these instruction assume you are using a Mac. "
},

{
    "location": "Installation/#Python-1",
    "page": "Installation",
    "title": "Python",
    "category": "section",
    "text": "Install python.  I recommend the anaconda distribution https://www.continuum.io/downloads.Once you install python, you are going to want two packages: a plotting package that Julia will use, and jupyter notebooks for interacting with Julia.  Install them likeconda install matplotlib\nconda install mathjax\nconda install jupyter"
},

{
    "location": "Installation/#Julia-1",
    "page": "Installation",
    "title": "Julia",
    "category": "section",
    "text": "You can get Julia from  http://julialang.org/.   If you are using a Mac, you may wish to create a symnolic link to the Julia executable so that you can call it from a terminal.  For example, you can do this like:cd /usr/local/bin/\nln -s julia /Applications/Julia-0.5.app/Contents/Resources/julia/bin/juliaOnce you have this, you will want Julia notebooks.  To install this, run julia and typejulia> Pkg.add(\"IJulia\")\njulia> using IJuliaThis will install the package, and put the current julia kernel into jupyter.  In the future, you can launch the Julia notebook by typing (in a terminal)jupyter notebookMost users will also want to install PyPlot, if you did not already. To do that, typePkg.add(\"PyPlot\")"
},

{
    "location": "Installation/#Laplacians-1",
    "page": "Installation",
    "title": "Laplacians",
    "category": "section",
    "text": "In theory, all you need to do now is type eitherjulia> Pkg.add(\"Laplacians\")To use the package, you then typejulia> using LaplaciansThe one catch is with the functions for drawing graphs.  These require PyPlot.  If you did not install it before typing Pkg.add(\"PyPlot\"), then you can either install it now or disable the plotting routines in Laplacians.If you do not want to load PyPlot, then either set the environment variable LAPLACIANS_NOPLOT to true in bash, like$ export LAPLACIANS_NOPLOT=trueor, set the variable inside Julia, like julia> LAPLACIANS_NOPLOT = truebefore typing using Laplacians. Similarly, you can avoid loading PyAmg by settingjulia> LAPLACIANS_NOAMG = true(note that these are not the same variable: the environment variable in Julia is available as ENV[\"LAPLACIANS_NOPLOT\"]. Actually, defining these variables to anything will have the same effect.  So, setting them to false has the same effect as setting them to true.To see if Laplacians is working, try typinga = chimera(100,6)\nspectral_drawing(a)ora = generalizedNecklace(grid2(6),grid2(3),2)\nspectral_drawing(a)"
},

{
    "location": "Examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "Examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "The following are links to html files of Julia notebooks. These notebooks are also in the notebook directory, and can be open there so that you can run the code live. You should be able to find them under ~/.julia/v0.4/Laplacians/notebooks.FirstNotebook\nSolvers\nSampler\nLocalClustering\nLocalClustering Statistics"
},

{
    "location": "CSCgraph/#",
    "page": "Sparse matrices as graphs",
    "title": "Sparse matrices as graphs",
    "category": "page",
    "text": ""
},

{
    "location": "CSCgraph/#Using-sparse-matrices-as-graphs-1",
    "page": "Sparse matrices as graphs",
    "title": "Using sparse matrices as graphs",
    "category": "section",
    "text": "The routines deg, nbri and weighti will let you treat a sparse matrix like a graph.deg(graph, u) is the degree of node u. nbri(graph, u, i) is the ith neighbor of node u. weighti(graph, u, i) is the weight of the edge to the ith neighbor of node u.Note that we start indexing from 1.For example, to iterate over the neighbors of node v,   and play with the attached nodes, you could write code like:  for i in 1:deg(mat, v)\n     nbr = nbri(mat, v, i)\n     wt = weighti(mat, v, i)\n     foo(v, nbr, wt)\n  endBut, this turns out to be much slower than working with the structure directly, like  for ind in mat.colptr[v]:(mat.colptr[v+1]-1)\n      nbr = mat.rowval[ind]\n      wt = mat.nzval[ind]\n      foo(v, nbr, wt)\n  end[ ] Maybe we can make a macro to replace those functions.  It could be faster and more readable."
},

{
    "location": "CSCgraph/#The-SparseMatrixCSC-data-structure-1",
    "page": "Sparse matrices as graphs",
    "title": "The SparseMatrixCSC data structure",
    "category": "section",
    "text": "You can explore what is going on with the data structure by looking at some examples.  For example, here is a randomly weighted complete graph on 4 vertices, first displayed as a matrix:gr = round(10*uniformWeight(completeGraph(4)))\n\n4x4 sparse matrix with 12 Float64 entries:\n	[2, 1]  =  3.0\n	[3, 1]  =  3.0\n	[4, 1]  =  6.0\n	[1, 2]  =  3.0\n	[3, 2]  =  1.0\n	[4, 2]  =  2.0\n	[1, 3]  =  3.0\n	[2, 3]  =  1.0\n	[4, 3]  =  7.0\n	[1, 4]  =  6.0\n	[2, 4]  =  2.0\n	[3, 4]  =  7.0\n	\nfull(gr)\n\n4x4 Array{Float64,2}:\n 0.0  3.0  3.0  6.0\n 3.0  0.0  1.0  2.0\n 3.0  1.0  0.0  7.0\n 6.0  2.0  7.0  0.0To see the underlying data structure, use fieldnames.fieldnames(gr)\n\n5-element Array{Symbol,1}:\n :m     \n :n     \n :colptr\n :rowval\n :nzval m and n are the dimensions of the matrix. The entries of the matrix are stored in nzval. colptr[i] is the index in nzval of the first nonzero entry in column i.  rowval tells you which rows in each column are nonzero. The indices of the nonzero entries in column i are stored in  rowval[colptr[i]] through rowval[colptr[i+1]-1].gr.colptr \n\n5-element Array{Int64,1}:\n  1\n  4\n  7\n 10\n 13\n \n [gr.rowval gr.nzval]\n \n 12x2 Array{Float64,2}:\n 2.0  3.0\n 3.0  3.0\n 4.0  6.0\n 1.0  3.0\n 3.0  1.0\n 4.0  2.0\n 1.0  3.0\n 2.0  1.0\n 4.0  7.0\n 1.0  6.0\n 2.0  2.0\n 3.0  7.0"
},

{
    "location": "usingSolvers/#",
    "page": "Solving Linear Equations",
    "title": "Solving Linear Equations",
    "category": "page",
    "text": "[TOC]"
},

{
    "location": "usingSolvers/#Solving-linear-equations-in-Laplacians-and-SDD-matrices-1",
    "page": "Solving Linear Equations",
    "title": "Solving linear equations in Laplacians and SDD matrices",
    "category": "section",
    "text": "The main purpose of this package is to experiment with the implementation of algorithms for solving systems of linear equations in Laplacian and symmetric, diagonally dominant, M-matrices (SDDM).At present, the fastest solver in this package for Laplacians is approxchol_lap. For SDDM systems, one should use approxchol_sddm.  Here is a quick demo.  Read more for other solvers and other options you can pass to the solvers.julia> a = grid3(50); # an adjacency matrix\njulia> la = lap(a); # it\'s Laplacian\njulia> sol = approxchol_lap(a); # a solver for la\njulia> b = randn(size(la,1)); b = b - mean(b); # a right-hand-side\njulia> x = sol(b); # the solution\njulia> norm(la*x-b) / norm(b)\n5.911931368666469e-7\njulia> x = sol(b, tol=1e-12); # a higher accuracy solution\njulia> norm(la*x-b) / norm(b)\n7.555529748070115e-11\njulia> x = sol(b, tol=1e-1, verbose=true); # faster, lower accuracy, with info\nPCG stopped after: 0.022 seconds and 3 iterations with relative error 0.07929402690389374.\n\njulia> sddm = copy(la); # doing it with a SDDM matrix\njulia> sddm[1,1] += 1;\njulia> sol = approxchol_sddm(sddm, verbose=true); # solver, with output\nUsing greedy degree ordering. Factorization time: 0.7143130302429199\nRatio of operator edges to original edges: 2.1120548223350255\nratio of max to min diagonal of laplacian : 6.0\nSolver build time: 0.747 seconds.\n\njulia> x = sol(b, verbose=false); # a solve, supressing output\njulia> norm(sddm*x - b) / norm(b)\n8.739618692868002e-7We recall that a matrix $ L $ is a Laplacian matrix if:It is symmetric,\nits off-diagonal entries are non-positive, and\nall of its row sums are 0.These conditions imply that the diagonal entries are non-negative, and that the matrix is singular.  So, we only hope to solve equations of the form  $ Lx = b $ when b is in the span of the matrix.  When the graph of the nonzero entries of the matrix is connected, this is precisely when the sum of the entries in $ b $ is zero.  Laplacian matrices are always positive semidefinite.A matrix $ M $ is a symmetric M-matrix if:It is symmetric,\nits off-diagonal entries are non-positive, and\nit is positive definite.A matrix symmetric $ M $ is diagonally dominant if each of its diagonals is at least the sum of the absolute values of the off-diagonal entries in its row.  A Laplacians matrix is diagonally dominant.  A diagonally dominant matrix is always positive semidefinite.A SDDM matrix (symmetric, diagonally-dominant M-matrix) is a matrix that is both diagonally dominant and an M-matrix.  You may think of a SDDM matrix as a Laplacian plus a non-negative, non-zero, diagonal matrix.  However, this is only guaranteed to produce a SDDM matrix when the graph underlying the Laplacian is connected.Laplacians.jl contains code for solving systems of linear equations in both Laplacian and SDDM matrices.  In fact, these problems are equivalent.  So, usually a solver for one type of system is implemented, and then wrapped to solve the other. The same ideas can be used to solve systems of equations in SDD matrices (the off-diagonals can be positive or negative), but a wrapper for these has not yet been written."
},

{
    "location": "usingSolvers/#Harmonic-Interpolation-1",
    "page": "Solving Linear Equations",
    "title": "Harmonic Interpolation",
    "category": "section",
    "text": "One of the main reasons to solve Laplacian and SDDM systems is to interpolate harmonic functions on graphs.  In an unweighted graph, these have the property that the value at every vertex is the average of the values of its neighbors.  To make this sensible, some values must be fixed.For example, below we fit a harmonic function on the 4-by-4 grid. We fix the values of vertices 1, 4, and 16 to 0.0, 0.5, and 2.0, respectively.  We then show the results by forcing them into a 4-by-4 grid.julia> a = grid2(4);\njulia> S = [1; 4; 16];\njulia> vals = [0; 0.5; 2];\njulia> x = harmonic_interp(a, S, vals);\njulia> reshape(x,4,4)\n4×4 Array{Float64,2}:\n 0.0       0.460942  0.749255  0.903101\n 0.398927  0.633572  0.883721  1.05695\n 0.563208  0.790698  1.09511   1.38402\n 0.5       0.8709    1.322     2.0     "
},

{
    "location": "usingSolvers/#The-Solver-Interface-1",
    "page": "Solving Linear Equations",
    "title": "The Solver Interface",
    "category": "section",
    "text": "All of the SDDM solvers take the SDDM matrix as input.All of the Laplacian solvers take the adjacency matrix of the underlying graph as input.To solve a system of linear equations, one first passes the matrix defining the system to a linear equation solving algorithm.  This will return a function that solves systems of linear equations in that matrix.  For example,julia> n = 1000;\njulia> a = wted_chimera(n);  # produces a graph, as a sparse adjacency matrix\njulia> b = randn(n);\njulia> b = b - mean(b); # so there is a solution\njulia> f = chol_lap(a)\n(::#79) (generic function with 1 method)\njulia> x = f(b);\njulia> la = lap(a);  # construct the Laplacian of a\njulia> norm(la*x-b)\n2.9565023548855584e-13All of the solvers take the following keyword arguments. This means that they are optional, and will be set to their default values if not specified.tol : the relative accuracy required: $ \\| M x - b \\| / \\| b \\| $.\nmaxits : the maximum number of iterations, for iterative algorithms.\nmaxtime : quit if it takes more than this many seconds.  Not all routines obey this, but they try.\nverbose : set to true to display diagnostic information.\npcgIts : If the algorithm is iterative, this allows it to return the number of iterations it performed.  If pcgIts is an array of positive length, then its first entry is set to the number of iterations.  Where verbose prints this information, pcgIts allows it to be returned to other code.  To disable this set pcgIts to a zero length array, like Int[].Most of the solvers are iterative, exploiting the preconditioned conjugate gradient. These are the solvers for which maxits, maxtime and pcgIts make the most sense.  Some solvers, like Cholesky factorization, just ignore these parameters.All of these parameters may be set in the call that constructs f.  They may then be over-ridden by again setting them in the call to f. Let\'s see how this works when using the conjugate gradient.julia> f = cgLapSolver(a, tol=1e-2, verbose=true)\n(::f) (generic function with 1 method)\njulia> x = f(b);\nCG BLAS stopped after: 78 iterations with relative error 0.009590493139133275.\njulia> norm(la*x-b)/norm(b)\n0.00959049313913375\n\njulia> pcgIts = [0]\n1-element Array{Int64,1}:\n 0\njulia> x = f(b,verbose=false, pcgIts=pcgIts);\njulia> pcgIts\n1-element Array{Int64,1}:\n 78\n\njulia> x = f(b,verbose=true, maxits=50);\nCG BLAS stopped after: 50 iterations with relative error 0.050483096216933886.\n\njulia> x = f(b, tol=1e-4);\nCG BLAS stopped after: 131 iterations with relative error 8.886882933346416e-5.\njulia> norm(la*x-b)/norm(b)\n8.886882933294668e-5For some experiments with solvers, including some of those below, look at the notebook Solvers.ipynb.In the following, we document many of the solvers that have been implemented in this package."
},

{
    "location": "usingSolvers/#Cholesky-Factorization-1",
    "page": "Solving Linear Equations",
    "title": "Cholesky Factorization",
    "category": "section",
    "text": "Cholesky factorization, the version of Gaussian Elimination for symmetric matrices, should be the first solver you try.  It will be very fast for matrices of dimension less than 1000, and for larger matrices coming from two-dimensional problems.You can compute a cholesky factor directly with cholfact.  It does  more than just compute the factor, and it saves its result in a data structure that implements \\.  It uses SuiteSparse by Tim Davis.Here is an example of how you would use it to solve a general non-singular linear system.a = grid2(5)\nla = lap(a)\nsddm = copy(la)\nsddm[1,1] = sddm[1,1] + 1\nF = cholfact(sddm)\n\nn = size(a)[1]\nb = randn(n)\nx = F \\ b\nnorm(sddm*x-b)\n\n 	1.0598778281116327e-14As cholfact does not satisfy our interface, we wrap it in a routine chol_sddm that does.To solve systems in Laplacian matrices, use chol_lap.  Recall that this routine should be passed the adjacency matrix.f = chol_lap(a)\nb = randn(n);\nb = b - mean(b);\nnorm(la*f(b) - b)\n	2.0971536951312585e-15"
},

{
    "location": "usingSolvers/#CG-and-PCG-1",
    "page": "Solving Linear Equations",
    "title": "CG and PCG",
    "category": "section",
    "text": "We have written implementations of Conjugate Gradient (CG) and Preconditioned Conjugate Gradient (PCG) that satisfy the interface. These routines use BLAS when possible, and slower routines when dealing with data types that BLAS cannot handle.  seed!(1)\nn = 50\nM = randn(n,n); M = M * M\';\nb = randn(n)\nx = cg(M,b,maxits=100,verbose=true);\nCG BLAS stopped after: 66 iterations with relative error 2.0166243927814765e-7.\n\nbbig = convert(Array{BigFloat,1},b)\nxbig = cg(M,bbig,maxits=100,tol=1e-30)\nCG Slow stopped after: 50 iterations with relative error 2.18672511297479336887519117065525148757254642683072581090418060286711737398731e-38.\n\nnorm(M*xbig - bbig)\n1.605742093628722039938504001423963138146137896744531914963345296279741402982296e-37To create a function f that uses cg to solve systems in M, use cgSolver.  For Laplacians, use cgLapSolver.julia> n = 1000;\njulia> a = wted_chimera(n,1);\njulia> f = Laplacians.cgLapSolver(a,maxits=100);\n\njulia> b = randn(n);\njulia> b = b - mean(b);\njulia> x = f(b,verbose=true);\nCG BLAS stopped after: 100 iterations with relative error 0.012102058751548373.\n\n\njulia> la = lap(a);\njulia> sddm = copy(la);\njulia> sddm = sddm + spdiagm(rand(n)/100);\njulia> g = cgSolver(sddm,verbose=true)\n(::f) (generic function with 1 method)\n\njulia> x = g(b);\nCG BLAS stopped after: 253 iterations with relative error 7.860172210007891e-7.PCG also takes as input a preconditioner.  This should be a function.  Here is an example of how one might construct and use a diagonal preonditioner.  To motivate this, I will use a grid with highly varying weights on edges.seed!(1)\na = mapweight(grid2(200),x->1/(rand(1)[1]));\nla = lap(a)\nn = size(la)[1]\nb = randn(n)\nb = b - mean(b);\n\nd = diag(la)\nprec(x) = x ./ d\n@time x = pcg(la,b,prec,maxtime=1,tol=1e-2,verbose=true);\n\nPCG BLAS stopped at maxtime.\nPCG BLAS stopped after: 530 iterations with relative error 0.07732478003311881.\n  1.007756 seconds (10.32 k allocations: 648.525 MB, 9.69% gc time)\n\n@time x = pcg(la,b,prec,maxtime=3,tol=1e-2,verbose=true);\nPCG BLAS stopped after: 1019 iterations with relative error 0.009984013184429813.\n  2.086828 seconds (19.57 k allocations: 1.216 GB, 9.92% gc time)Without the preconditioner, CG takes much longer on this example.@time x = cg(la,b,tol=1e-2,maxtime=10,verbose=true);\n\nCG BLAS stopped at maxtime.\nCG BLAS stopped after: 8879 iterations with relative error 0.054355534834831624.\n 10.001998 seconds (97.91 k allocations: 2.649 GB, 4.48% gc time)pcgSolver creates a function that uses the preconditioner to solve systems in the matrix.f = pcgSolver(la,prec)\n@time x = f(b,maxtime=3,tol=1e-2,verbose=true);\nPCG BLAS stopped after: 1019 iterations with relative error 0.009984013184429813.\n  1.892217 seconds (19.58 k allocations: 1.216 GB, 9.47% gc time)pcgLapSolver uses the Laplacian of one matrix as a preconditioner for the first.  It solves systems of linear equations in the preconditioner by Cholesky factorization.  It performs the Cholesky factorization when pcgLapSolver is called.  This is why we do the work of creating f only once.  Here is an example using a Low-Stretch Spanning Tree preconditioner.@time t = akpw(a)\n  0.210467 seconds (1.43 M allocations: 91.226 MB, 19.23% gc time)\n\n@time f = pcgLapSolver(a,t)\n  0.160210 seconds (288 allocations: 28.076 MB, 72.28% gc time)\n\n@time x = f(b,maxtime=3,tol=1e-2,verbose=true);\nPCG BLAS stopped after: 260 iterations with relative error 0.009864463201800925.\n  1.014897 seconds (28.02 k allocations: 1.008 GB, 9.81% gc time)"
},

{
    "location": "usingSolvers/#Low-Stretch-Spanning-Trees-1",
    "page": "Solving Linear Equations",
    "title": "Low-Stretch Spanning Trees",
    "category": "section",
    "text": "In order to make preconditioners, we will want low-stretch spanning trees.  We do not yet have any code that is guaranteed to produce these.  Instead, we supply three heuristics: akpw which is inspired by the algorith of Alon, Karp, Peleg and West, and  randomized versions of Prim and Kruskal\'s algorithm. randishKruskal samples the remaining edges with probability proportional to their weight.  randishPrim samples edges on the boundary while using the same rule.  We recommend using akpw.See Low Stretch Spanning Trees to learn more about these."
},

{
    "location": "usingSolvers/#Augmented-Spanning-Tree-Preconditioners-1",
    "page": "Solving Linear Equations",
    "title": "Augmented Spanning Tree Preconditioners",
    "category": "section",
    "text": "These are obtained by constructing a spanning tree of a graph, and then adding back some more edges from the graph.  The tree should have low stretch.  The edges to add back are chosen at random with probabilities proportional to their stretches.These are implemented in the routinesaugTreeSddm, for SDDM matrices\naugTreeLap\naugTreePrecon\naugmentTree"
},

{
    "location": "usingSolvers/#The-solvers-of-Koutis,-Miller-and-Peng.-1",
    "page": "Solving Linear Equations",
    "title": "The solvers of Koutis, Miller and Peng.",
    "category": "section",
    "text": "Solvers inspired by the algorithm from \"Approaching optimality for solving SDD systems\" by Koutis, Miller, and Peng, <i>SIAM Journal on Computing</i>, 2014.KMPLapSolver\nKMPSDDMSolver"
},

{
    "location": "usingSolvers/#Sampling-Solvers-of-Kyng-and-Sachdeva-1",
    "page": "Solving Linear Equations",
    "title": "Sampling Solvers of Kyng and Sachdeva",
    "category": "section",
    "text": "These are inspired by the paper \"Approximate Gaussian Elimination for Laplacians: Fast, Sparse, and Simple\" by Rasmus Kyng and Sushant Sachdeva, FOCS 2016.These first two follow that paper reasonably closely.samplingSDDMSolver\nsamplingLapSolverThe following is a modification of the algorithm that eliminates edges one at a time.  The code is by Daniel Spielman.  The algorithm has not yet been analyzed.  It is presently the fastest in this package.approxchol_lap"
},

{
    "location": "usingSolvers/#Algebraic-Multigrid-1",
    "page": "Solving Linear Equations",
    "title": "Algebraic Multigrid",
    "category": "section",
    "text": "This is an interface to the algebraic multigrid solvers from the PyAMG package.AMGLapSolver\nAMGSolver, for SDDM systems."
},

{
    "location": "usingSolvers/#Solvers-from-Matlab-1",
    "page": "Solving Linear Equations",
    "title": "Solvers from Matlab",
    "category": "section",
    "text": "The MATLAB.jl package allows Julia to call routines from Matlab, provided you have Matlab installed.  It does this in a very efficient fashion: it starts up the Matlab process when you type using MATLAB, and then communicates with it.  So, we have wrapped some solvers from Matlab so that they obey the same interface.These are not part of the Laplacians module, but are included in the package under src/matlabSolvers.jl.  To include them, typeinclude(string(Pkg.dir(\"Laplacians\") , \"/src/matlabSolvers.jl\"))We provide the docstrings for these here."
},

{
    "location": "usingSolvers/#Incomplete-Cholesky-Factorizations-1",
    "page": "Solving Linear Equations",
    "title": "Incomplete Cholesky Factorizations",
    "category": "section",
    "text": "These use the no-fill incomplete Cholesky factorizations implemented in Matlab.  They first order the vertices by the symrcm ordering.The solvers are:f = matlab_ichol_sddm(sddm; tol, maxtime, maxits, pctIts, verbose)\nf = matlab_ichol_lap(A; tol, maxtime, maxits, pctIts, verbose)A routine that just wraps the function that solves equations in the preconditioner is provided as well:f = matlab_ichol(sddm)"
},

{
    "location": "usingSolvers/#Koutis\'s-Combinatorial-Multigrid-(CMG)-1",
    "page": "Solving Linear Equations",
    "title": "Koutis\'s Combinatorial Multigrid (CMG)",
    "category": "section",
    "text": "You must have installed Yiannis Koutis\'s Combinatorial Multigrid Code, and it must be on Matlab\'s default path.  As this code returns a function rather than a preconditioner, it would be inefficient to make it use our PCG code and satisfy our interface.  So, it does not.x = matlabCmgSolver(mat, b; tol::Real=1e-6, maxits=10000)The matrix mat can either be SDDM or a Laplacian.  This solves the system in b.If you need to specify the solver separately from b, you can callx = matlabCmgSolver(mat; tol::Real=1e-6, maxits=10000)or, for the Laplacians of the adjacency matrix A,x = matlabCmgLap(A; tol::Real=1e-6, maxits=10000)However, this does not create the solver.  It merely returns a call to the previous routine."
},

{
    "location": "LSST/#",
    "page": "Low Stretch Spanning Trees",
    "title": "Low Stretch Spanning Trees",
    "category": "page",
    "text": ""
},

{
    "location": "LSST/#Low-Stretch-Spanning-Trees-1",
    "page": "Low Stretch Spanning Trees",
    "title": "Low Stretch Spanning Trees",
    "category": "section",
    "text": "We have implemented a variant of the algorithm of Alon, Karp, Peleg and West for computing low stretch spanning trees.  It is called akpw.  For unweighted graphs we provide a faster variant called akpwU.  If you require faster algorithms, with possibly higher average stretch, you can try the heuristics randishPrim or randishKruskal.You can compute the stretch of a graph with respect to a spanning tree with the routine comp_stretches.  It returns a sparse matrix with one entry for each edge in the graph giving its stretch.  For example:julia> graph = grid2(4)\n16x16 sparse matrix with 48 Float64 entries:\n\njulia> tree = akpwU(graph)\n16x16 sparse matrix with 30 Float64 entries:\n\njulia> st = comp_stretches(tree,graph)\n16x16 sparse matrix with 48 Float64 entries:\n	[2 ,  1]  =  1.0\n	[5 ,  1]  =  1.0\n	[1 ,  2]  =  1.0\n	[3 ,  2]  =  1.0\n	[6 ,  2]  =  3.0\n	[2 ,  3]  =  1.0\n	[4 ,  3]  =  1.0\n	[7 ,  3]  =  1.0\n	[3 ,  4]  =  1.0\n	[8 ,  4]  =  3.0\n	[1 ,  5]  =  1.0\n	[6 ,  5]  =  5.0\n	⋮\n	[8 , 12]  =  3.0\n	[11, 12]  =  1.0\n	[16, 12]  =  1.0\n	[9 , 13]  =  1.0\n	[14, 13]  =  3.0\n	[10, 14]  =  1.0\n	[13, 14]  =  3.0\n	[15, 14]  =  3.0\n	[11, 15]  =  1.0\n	[14, 15]  =  3.0\n	[16, 15]  =  3.0\n	[12, 16]  =  1.0\n	[15, 16]  =  3.0Here is an example demonstrating the average stretches and times taken by these algorithms on a large graph.\njulia> graph = chimera(1000000,1);\n\njulia> @time tree = akpw(graph);\n  5.700285 seconds (16.10 M allocations: 1.263 GB, 11.16% gc time)\n\njulia> aveStretch = sum(comp_stretches(tree,graph))/nnz(graph)\n8.793236275779616\n\njulia> @time tree = randishPrim(graph);\n  3.736225 seconds (3.21 M allocations: 566.887 MB, 6.40% gc time)\n\njulia> aveStretch = sum(comp_stretches(tree,graph))/nnz(graph)\n10.800094649795756\n\njulia> @time tree = randishKruskal(graph);\n  2.515443 seconds (3.21 M allocations: 423.529 MB, 4.35% gc time)\n\njulia> aveStretch = sum(comp_stretches(tree,graph))/nnz(graph)\n37.819948689847564\nOf course, you can get very different results on very different graphs.  But, the ordering of these algorithms respect to time and stretch will usually follow this pattern."
},

{
    "location": "Developing/#",
    "page": "Developing",
    "title": "Developing",
    "category": "page",
    "text": ""
},

{
    "location": "Developing/#Developing-Laplacians.jl-1",
    "page": "Developing",
    "title": "Developing Laplacians.jl",
    "category": "section",
    "text": ""
},

{
    "location": "Developing/#Learn-to-use-git-1",
    "page": "Developing",
    "title": "Learn to use git",
    "category": "section",
    "text": "If you don\'t know anything about git, then just know that you should make a branch for you own code.  Typegit checkout -b MyNameMake sure that your .gitignore file contains the lines*~\n*#\n.ipynb_*\n.DS_Store\n*.cov\ndocs/build\ndocs/siteNow, read about Git.  I recommend the book Pro Git, which is available online for free.\nStop thinking about Git like subversion or dropbox.\nThe master branch will be the one for public consumption. It should (mostly) work.\nYou should also read the section of the Julia docs about building packages."
},

{
    "location": "Developing/#Tests-1",
    "page": "Developing",
    "title": "Tests",
    "category": "section",
    "text": "Every piece of code should have a test in the \"tests\" directory.  Ideally, it should have enough tests to run every line of the code.  To run all the tests from the command prompt, typejulia -e \'Pkg.test(\"Laplacians\")\'"
},

{
    "location": "Developing/#Fast-code?-1",
    "page": "Developing",
    "title": "Fast code?",
    "category": "section",
    "text": "Just go for it. Don\'t worry about writing fast code at first. Just get it to work. We can speed it up later.Within some of the files, I am keeping old, unoptimized versions of code around for comparison (and for satisfaction).  I will give them the name \"XSlow\""
},

{
    "location": "Developing/#Documentation-1",
    "page": "Developing",
    "title": "Documentation",
    "category": "section",
    "text": "This documentation is still very rough. It is generated by a combination of Markdown and semi-automatic generation, using the Documenter.jl package.  The one annoying feature of this package is that it will not allow the inclusion of a docstring on more than one page.  I don\'t know why.The steps to generate and improve the documentation are:Edit Markdown files in the docs directory.  For example, you could use MacDown to do this.\nIf you want to add a new page to the documention, create one.  Edit the file mkdocs.yml so show where it should appear.\nAdd docstrings to everything that needs it, and in particular to the routines you create.  The API is built from the docstrings. \nRun julia make.jl; mkdocs build in the docs directory to generate the documentation from the Markdown.  This will generate a local copy of the documentation that you can use for reference.\nWARNING: You should not include any pages that are generated in the git repository.  So, make sure that your .gitignore file contains the line docs/build and docs/site.When you push to the master branch on GitHub, Travis will automatically build and update the docs.  <b>DO NOT RUN mkdocs gh-deploy</b>If you create a Julia notebook that you would like to include as documentation.   You should  put it in the notebooks directory (.julia/v0.5/Laplacians/notebooks) and then link to it\'s page on GitHub.  While it seems that one should convert it to html (and one can), and then include it in MkDocs, MkDocs does something funny to the resulting html that does not look nice."
},

{
    "location": "Developing/#Parametric-Types-1",
    "page": "Developing",
    "title": "Parametric Types",
    "category": "section",
    "text": "A sparse matrix has two types associated with it: the types of its indices (some sort of integer) and the types of its values (some sort of number).  Most of the code has been written so that once these types are fixed, the type of everything else in the function has been too.  This is accomplished by putting curly braces after a function name, with the names of the types that we want to use in the braces.  For example,shortestPaths{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti}, start::Ti)Tv, sometimes written Tval denotes the types of the values, and Ti or Tind denotes the types of the indices.  This function will only be called if the node from which we compute the shortest paths, start is of type Ti.  Inside the code, whenever we write something like pArray = zeros(Ti,n), it creates an array of zeros of type Ti.  Using these parameteric types is much faster than leaving the types unfixed."
},

{
    "location": "Developing/#Data-structures:-1",
    "page": "Developing",
    "title": "Data structures:",
    "category": "section",
    "text": "IntHeap a heap that stores small integers (like indices of nodes in a graph) and that makes deletion fast.  Was much faster than using Julia\'s more general heap."
},

{
    "location": "Developing/#Interface-issue:-1",
    "page": "Developing",
    "title": "Interface issue:",
    "category": "section",
    "text": "There are many different sorts of things that our code could be passing around.  For example, kruskal returns a graph as a sparse matrix.  But, we could use a format that is more specialized for trees, like the RootedTree type.  At some point, when we optimize code, we will need to figure out the right interfaces between routines.  For example, some routines symmetrize at the end.  This is slow, and should be skipped if not necessary.  It also doubles storage."
},

{
    "location": "Developing/#Integration-with-other-packages.-1",
    "page": "Developing",
    "title": "Integration with other packages.",
    "category": "section",
    "text": "There are other graph packages that we might want to sometimes use.Graphs.jl : I found this one to be too slow and awkward to be useful.\nLightGraphs.jl : this looks more promising.  We will have to check it out.  It is reasonably fast, and the code looks pretty."
},

{
    "location": "graphGenerators/#",
    "page": "generators",
    "title": "generators",
    "category": "page",
    "text": ""
},

{
    "location": "graphGenerators/#Generators-1",
    "page": "generators",
    "title": "Generators",
    "category": "section",
    "text": "Laplacians.jl implements generators for many standard graphs. The chimera and wted_chimera generators are designed to stress code by combining these standard graphs in tricky ways.  While no one of these graphs need be a hard case for any application, the goal is for these generators to explore the space of graphs in such a way that running on many of them should exercise your code.chimera(n) generates a random chimera graph. chimera(n,k) first sets the seed of the psrg to k. In this way, it generates the kth chimera graph, and messes with your psrg. wted_chimera is similar, but it generates weighted graphs."
},

{
    "location": "graphGenerators/#Laplacians.ErdosRenyi-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.ErdosRenyi",
    "category": "method",
    "text": "graph = ErdosRenyi(n::Integer, m::Integer; ver=Vcur)\n\nGenerate a random graph on n vertices with m edges. The actual number of edges will probably be smaller, as we sample with replacement\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.ErdosRenyiCluster-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.ErdosRenyiCluster",
    "category": "method",
    "text": "graph = ErdosRenyiCluster(n::Integer, k::Integer; ver=Vcur)\n\nGenerate an ER graph with average degree k, and then return the largest component. Will probably have fewer than n vertices. If you want to add a tree to bring it back to n, try ErdosRenyiClusterFix.\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.ErdosRenyiClusterFix-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.ErdosRenyiClusterFix",
    "category": "method",
    "text": "graph = ErdosRenyiClusterFix(n::Integer, k::Integer; ver=Vcur)\n\nLike an Erdos-Renyi cluster, but add back a tree so it has n vertices\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.chimera-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.chimera",
    "category": "method",
    "text": "graph = chimera(n::Integer, k::Integer; verbose=false, ver=Vcur)\n\nBuilds the kth chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.chimera-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.chimera",
    "category": "method",
    "text": "graph = chimera(n::Integer; verbose=false, ver=Vcur)\n\nBuilds a chimeric graph on n vertices. The components come from pureRandomGraph, connected by joinGraphs, productGraph and generalizedNecklace\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.complete_binary_tree-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.complete_binary_tree",
    "category": "method",
    "text": "graph = completeBinaryTree(n::Int64)\n\nThe complete binary tree on n vertices\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.complete_bipartite_graph-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.complete_bipartite_graph",
    "category": "method",
    "text": "graph = complete_bipartite_graph(n)\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.complete_graph-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.complete_graph",
    "category": "method",
    "text": "graph = complete_graph(n)\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.empty_graph-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.empty_graph",
    "category": "method",
    "text": "ijv = empty_graph(n)\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.generalized_ring-Union{Tuple{T}, Tuple{T,Array{T,N} where N}} where T<:Integer",
    "page": "generators",
    "title": "Laplacians.generalized_ring",
    "category": "method",
    "text": "graph = generalized_ring(n, gens)\n\nA generalization of a ring graph. The vertices are integers modulo n. Two are connected if their difference is in gens. For example,\n\ngeneralized_ring(17, [1 5])\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.grid2-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.grid2",
    "category": "method",
    "text": "graph = grid2(n::Int64, m::Int64; isotropy=1)\n\nAn n-by-m grid graph.  iostropy is the weighting on edges in one direction.\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.grid2coords-Tuple{Int64,Int64}",
    "page": "generators",
    "title": "Laplacians.grid2coords",
    "category": "method",
    "text": "graph = grid2coords(n::Int64, m::Int64)\ngraph = grid2coords(n::Int64)\n\nCoordinates for plotting the vertices of the n-by-m grid graph\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.grid3-Tuple{Integer,Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.grid3",
    "category": "method",
    "text": "graph = grid3(n1, n2, n3)\ngraph = grid3(n)\n\nAn n1-by-n2-by-n3 grid graph.\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.grown_graph-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.grown_graph",
    "category": "method",
    "text": "graph = grown_graph(n, k; ver=Vcur)\n\nCreate a graph on n vertices. For each vertex, give it k edges to randomly chosen prior vertices. This is a variety of a preferential attachment graph.\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.grown_graph_d-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.grown_graph_d",
    "category": "method",
    "text": "graph = grown_graph_d(n::Integer, k::Integer; ver=Vcur)\n\nLike a grownGraph, but it forces the edges to all be distinct. It starts out with a k+1 clique on the first k vertices\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.hypercube-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.hypercube",
    "category": "method",
    "text": "graph = hyperCube(d::Int64)\n\nThe d dimensional hypercube.  Has 2^d vertices and d*2^(d-1) edges.\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.path_graph-Tuple{Any}",
    "page": "generators",
    "title": "Laplacians.path_graph",
    "category": "method",
    "text": "graph = path_graph(n)\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.pref_attach-Tuple{Integer,Integer,Float64}",
    "page": "generators",
    "title": "Laplacians.pref_attach",
    "category": "method",
    "text": "graph = pref_attach(n::Int64, k::Int64, p::Float64; ver=Vcur)\n\nA preferential attachment graph in which each vertex has k edges to those that come before.  These are chosen with probability p to be from a random vertex, and with probability 1-p to come from the endpoint of a random edge. It begins with a k-clique on the first k+1 vertices.\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.pure_random_graph-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.pure_random_graph",
    "category": "method",
    "text": "graph = pure_random_graph(n::Integer; verbose=false, ver=Vcur)\n\nGenerate a random graph with n vertices from one of our natural distributions\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.rand_gen_ring-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.rand_gen_ring",
    "category": "method",
    "text": "graph = rand_gen_ring(n, k; verbose = false, ver=Vcur)\n\nA random generalized ring graph of degree k. Gens always contains 1, and the other k-1 edge types are chosen from an exponential distribution\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.rand_matching-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.rand_matching",
    "category": "method",
    "text": "graph = rand_matching(n::Integer; ver=Vcur)\n\nA random matching on n vertices\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.rand_regular-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.rand_regular",
    "category": "method",
    "text": "graph = rand_regular(n, k; ver=Vcur)\n\nA sum of k random matchings on n vertices\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.rand_weight-Tuple{Any}",
    "page": "generators",
    "title": "Laplacians.rand_weight",
    "category": "method",
    "text": "graph = randWeight(graph; ver=Vcur)\n\nApplies one of a number of random weighting schemes to the edges of the graph\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.ring_graph-Tuple{Any}",
    "page": "generators",
    "title": "Laplacians.ring_graph",
    "category": "method",
    "text": "graph = ring_graph(n)\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.semiwted_chimera-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.semiwted_chimera",
    "category": "method",
    "text": "graph = semiwted_chimera(n::Integer; verbose=false, ver=Vcur)\n\nA Chimera graph with some weights.  The weights just appear when graphs are combined. For more interesting weights, use wted_chimera\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.star_graph-Tuple{Any}",
    "page": "generators",
    "title": "Laplacians.star_graph",
    "category": "method",
    "text": "graph = star_graph(n)\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.wted_chimera-Tuple{Integer,Integer}",
    "page": "generators",
    "title": "Laplacians.wted_chimera",
    "category": "method",
    "text": "graph = wted_chimera(n::Integer, k::Integer; verbose=false, ver=Vcur)\n\nBuilds the kth wted chimeric graph on n vertices. It does this by resetting the random number generator seed. It should captute the state of the generator before that and then return it, but it does not yet.\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Laplacians.wted_chimera-Tuple{Integer}",
    "page": "generators",
    "title": "Laplacians.wted_chimera",
    "category": "method",
    "text": "graph = wted_chimera(n::Integer)\n\nGenerate a chimera, and then apply a random weighting scheme\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Random.randperm-Tuple{AbstractArray{T,2} where T}",
    "page": "generators",
    "title": "Random.randperm",
    "category": "method",
    "text": "graph = randperm(mat::AbstractMatrix)\n        randperm(f::Expr)\n\nRandomly permutes the vertex indices\n\n\n\n\n\n"
},

{
    "location": "graphGenerators/#Function-list-1",
    "page": "generators",
    "title": "Function list",
    "category": "section",
    "text": "Order = [:type, :function]\nPages   = [\"graphGenerators.md\"]Modules = [Laplacians]\nPages   = [\"graphGenerators.jl\"]\nPrivate = false"
},

{
    "location": "operators/#",
    "page": "operators",
    "title": "operators",
    "category": "page",
    "text": ""
},

{
    "location": "operators/#Operators-1",
    "page": "operators",
    "title": "Operators",
    "category": "section",
    "text": "Operators transform graphs to produce new graphs."
},

{
    "location": "operators/#Laplacians.adj-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "operators",
    "title": "Laplacians.adj",
    "category": "method",
    "text": "a,d = adj(sddm)\n\nCreate an adjacency matrix and a diagonal vector from an SDD M-matrix. That is, from a Laplacian with added diagonal weights\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.diagmat-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "operators",
    "title": "Laplacians.diagmat",
    "category": "method",
    "text": "d = diagmat(a)\n\nReturns the diagonal weighted degree matrix(as a sparse matrix) of a graph\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.disjoin-Tuple{SparseArrays.SparseMatrixCSC,SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.disjoin",
    "category": "method",
    "text": "graph = disjoin(a,b)\n\nCreate a disjoint union of graphs a and b,   with no edges between them.\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.edgeVertexMat-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.edgeVertexMat",
    "category": "method",
    "text": "U = edgeVertexMat(a)\n\nThe signed edge-vertex adjacency matrix\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.floatGraph-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.floatGraph",
    "category": "method",
    "text": "graph = floatGraph(a::SparseMatrixCSC)\n\nConvert the nonzero entries in a graph to Float64.\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.join_graphs!-Tuple{Laplacians.IJV,Laplacians.IJV,Integer}",
    "page": "operators",
    "title": "Laplacians.join_graphs!",
    "category": "method",
    "text": "graph = join_graphs!(a::IJV, b::IJV, k::Integer)\n\nCreate a disjoint union of graphs a and b,  and then put k random edges between them, merging b into a.\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.join_graphs-Union{Tuple{Tind}, Tuple{Tval}, Tuple{SparseArrays.SparseMatrixCSC{Tval,Tind},SparseArrays.SparseMatrixCSC{Tval,Tind},Integer}} where Tind where Tval",
    "page": "operators",
    "title": "Laplacians.join_graphs",
    "category": "method",
    "text": "graph = joinGraphs(a, b, k::Integer)\n\nCreate a disjoint union of graphs a and b,  and then put k random edges between them\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.lap-Tuple{Any}",
    "page": "operators",
    "title": "Laplacians.lap",
    "category": "method",
    "text": "l = lap(a)\n\nCreate a Laplacian matrix from an adjacency matrix. We might want to do this differently, say by enforcing symmetry\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.line_graph-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.line_graph",
    "category": "method",
    "text": "H = line_graph(G::SparseMatrixCSC)\n\nLet G = (V, E) be a graph. The line graph of G is the graph whose vertices are the edges of G in which two are connected if they share an endpoint in G. That is, (u, v),(w, z) is an edge of the line graph if one of {u, v} is the same as one of {w, z}\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.mapweight-Union{Tuple{Tind}, Tuple{Tval}, Tuple{SparseArrays.SparseMatrixCSC{Tval,Tind},Any}} where Tind where Tval",
    "page": "operators",
    "title": "Laplacians.mapweight",
    "category": "method",
    "text": "b = mapweight(a, x->rand())\n\nCreate a new graph that is the same as the original, but with f applied to each nonzero entry of a. For example, to make the weight of every edge uniform in [0,1], we could write.\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.plot_graph-NTuple{4,Any}",
    "page": "operators",
    "title": "Laplacians.plot_graph",
    "category": "method",
    "text": "plot_graph(gr,x,y,z;color=[0,0,1],dots=true,setaxis=true,number=false)\n\nPlots graph gr with coordinates (x,y,z)\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.plot_graph-Tuple{Any,Any,Any}",
    "page": "operators",
    "title": "Laplacians.plot_graph",
    "category": "method",
    "text": "plot_graph(gr,x,y;color=[0,0,1],dots=true,setaxis=true,number=false)\n\nPlots graph gr with coordinates (x,y)\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.power-Tuple{SparseArrays.SparseMatrixCSC,Int64}",
    "page": "operators",
    "title": "Laplacians.power",
    "category": "method",
    "text": "ap = power(a::SparseMatrixCSC, k::Int)\n\nReturns the kth power of a.\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.product_graph-Tuple{SparseArrays.SparseMatrixCSC,SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.product_graph",
    "category": "method",
    "text": "aprod = productGraph(a0, a1)\n\nThe Cartesian product of two graphs.  When applied to two paths, it gives a grid.\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.shortIntGraph-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.shortIntGraph",
    "category": "method",
    "text": "graph = shortIntGraph(a::SparseMatrixCSC)\n\nConvert the indices in a graph to 32-bit ints. This takes less storage, but does not speed up much.\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.spectral_coords-Tuple{Any}",
    "page": "operators",
    "title": "Laplacians.spectral_coords",
    "category": "method",
    "text": "x, y = spectral_coords(a)\n\nComputes the spectral coordinates of a graph. If more than 2 coords are desired, you can use\n\n    x, y, z = spectral_coords(a; k = 3)\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.spectral_drawing-Tuple{Any}",
    "page": "operators",
    "title": "Laplacians.spectral_drawing",
    "category": "method",
    "text": "spectral_drawing(a)\n\nComputes spectral coordinates, and then uses plot_graph to draw\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.subsampleEdges-Tuple{SparseArrays.SparseMatrixCSC,Float64}",
    "page": "operators",
    "title": "Laplacians.subsampleEdges",
    "category": "method",
    "text": "graph = subsampleEdges(a::SparseMatrixCSC, p::Float64)\n\nCreate a new graph from the old, but keeping edge edge with probability p\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.thicken-Tuple{SparseArrays.SparseMatrixCSC,Any}",
    "page": "operators",
    "title": "Laplacians.thicken",
    "category": "method",
    "text": "a_new = thicken(A,k)\n\nCreate a new graph with at least k times as many edges as A By connecting nodes with common neighbors at random. When this stops working (not enough new edges), repeat on the most recently produced graph. If k is too big, it is decreased so the average degree will not be pushed much above n/2.\n\nWhen called without k, it just runs thicken_once.\n\nFor example:\n\na = grid2(5)\na2 = thicken(a,3)\n(x,y) = grid2coords(5,5);\nplotGraph(a2,x,y)\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.thicken_once-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.thicken_once",
    "category": "method",
    "text": "a_new = thicken_once(a)\n\nCreates one edge for every vertex in a of degree > 1 by connecting two of its random neighbors. To use this to thicken a, return unweight(a + a_new).\n\na = grid2(5)\na2 = unweight(a + thicken_once(a))\n(x,y) = grid2coords(5,5);\nplotGraph(a2,x,y)\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.two_lift-Tuple{SparseArrays.SparseMatrixCSC,AbstractArray{Bool,1}}",
    "page": "operators",
    "title": "Laplacians.two_lift",
    "category": "method",
    "text": "graph = two_lift(a, flip::AbstractArray{Bool,1})\ngraph = two_lift(a)\ngraph = two_lift(a, k::Integer)\n\nCreats a 2-lift of a.  flip is a boolean indicating which edges cross. In the third version, k is the number of edges that cross.\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.uniformWeight!-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.uniformWeight!",
    "category": "method",
    "text": "uniformWeight!(a)\n\nSet the weight of every edge to random uniform [0,1]\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.unweight!-Union{Tuple{SparseArrays.SparseMatrixCSC{Tval,Tind}}, Tuple{Tind}, Tuple{Tval}} where Tind where Tval",
    "page": "operators",
    "title": "Laplacians.unweight!",
    "category": "method",
    "text": "unweight!(a)\n\nChange the weight of every edge in a to 1\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.unweight-Union{Tuple{SparseArrays.SparseMatrixCSC{Tval,Tind}}, Tuple{Tind}, Tuple{Tval}} where Tind where Tval",
    "page": "operators",
    "title": "Laplacians.unweight",
    "category": "method",
    "text": "wt1 = unweight(a)\n\nCreate a new graph in that is the same as the original, but with all edge weights 1\n\n\n\n\n\n"
},

{
    "location": "operators/#Laplacians.wtedEdgeVertexMat-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "operators",
    "title": "Laplacians.wtedEdgeVertexMat",
    "category": "method",
    "text": "U = wtedEdgeVertexMat(a)\n\nThe signed and weighted edge-vertex adjacency matrix, so U\'*U = L\n\n\n\n\n\n"
},

{
    "location": "operators/#Function-list-1",
    "page": "operators",
    "title": "Function list",
    "category": "section",
    "text": "Order = [:type, :function]\nPages   = [\"graphOps.md\"]Modules = [Laplacians]\nPages   = [\"graphOps.jl\"]\nPrivate = false"
},

{
    "location": "graphUtils/#",
    "page": "graphUtils",
    "title": "graphUtils",
    "category": "page",
    "text": ""
},

{
    "location": "graphUtils/#Laplacians.backIndices-Union{Tuple{Array{Array{Tuple{Tv1,Tv2},1},1}}, Tuple{Tv2}, Tuple{Tv1}} where Tv2 where Tv1",
    "page": "graphUtils",
    "title": "Laplacians.backIndices",
    "category": "method",
    "text": "Same as the above, but now the graph is in adjacency list form \n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Laplacians.backIndices-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "graphUtils",
    "title": "Laplacians.backIndices",
    "category": "method",
    "text": "Computes the back indices in a graph in O(M+N). works if for every edge (u,v), (v,u) is also in the graph \n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Laplacians.compConductance-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1}}} where Ti where Tv",
    "page": "graphUtils",
    "title": "Laplacians.compConductance",
    "category": "method",
    "text": "Returns the quality of the cut for a given graph and a given cut set s.   the result will be |outgoing edges| / min(|vertices in set|, |N - vertices in set|)\n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Laplacians.findEntries-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "graphUtils",
    "title": "Laplacians.findEntries",
    "category": "method",
    "text": "Similar to findnz, but also returns 0 entries that have an edge in the sparse matrix \n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Laplacians.flipIndex-Union{Tuple{SparseArrays.SparseMatrixCSC{Tval,Tind}}, Tuple{Tind}, Tuple{Tval}} where Tind where Tval",
    "page": "graphUtils",
    "title": "Laplacians.flipIndex",
    "category": "method",
    "text": "For a symmetric matrix, this gives the correspondance between pairs of entries in an ijv. So, ai[ind] = aj[flip[ind]].  For example, \n\n(ai,aj,av) = findnz(a);\nfl = flipIndex(a)\nind = 10\n@show backind = fl[10]\n@show [ai[ind], aj[ind], av[ind]]\n@show [ai[backind], aj[backind], av[backind]];\n\nbackind = fl[10] = 4\n[ai[ind],aj[ind],av[ind]] = [2.0,4.0,0.7]\n[ai[backind],aj[backind],av[backind]] = [4.0,2.0,0.7]\n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Laplacians.getObound-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1}}} where Ti where Tv",
    "page": "graphUtils",
    "title": "Laplacians.getObound",
    "category": "method",
    "text": "Computes the number of edges leaving s \n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Laplacians.getVolume-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1}}} where Ti where Tv",
    "page": "graphUtils",
    "title": "Laplacians.getVolume",
    "category": "method",
    "text": "Computes the volume of subset s in an unweighted graph G \n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Laplacians.setValue-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Ti,Ti,Tv}} where Ti where Tv",
    "page": "graphUtils",
    "title": "Laplacians.setValue",
    "category": "method",
    "text": "Sets the value of a certain edge in a sparse graph; value can be 0 without the edges dissapearing \n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Laplacians.wdeg-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Ti}} where Ti where Tv",
    "page": "graphUtils",
    "title": "Laplacians.wdeg",
    "category": "method",
    "text": "Finds the weighted degree of a vertex in the graph \n\n\n\n\n\n"
},

{
    "location": "graphUtils/#Graph-Utilities-1",
    "page": "graphUtils",
    "title": "Graph Utilities",
    "category": "section",
    "text": "These are utilities to facilitate the use of sparse matrices as graphs.Order = [:type, :function]\nPages   = [\"graphUtils.md\"]Modules = [Laplacians]\nPages   = [\"graphUtils.jl\"]\nPrivate = false"
},

{
    "location": "graphAlgs/#",
    "page": "graphAlgs",
    "title": "graphAlgs",
    "category": "page",
    "text": ""
},

{
    "location": "graphAlgs/#Graph-Algorithms-1",
    "page": "graphAlgs",
    "title": "Graph Algorithms",
    "category": "section",
    "text": "These are basic graph algorithms."
},

{
    "location": "graphAlgs/#Laplacians.biggestComp-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "graphAlgs",
    "title": "Laplacians.biggestComp",
    "category": "method",
    "text": "Return the biggest component in a graph, as a graph\n\n\n\n\n\n"
},

{
    "location": "graphAlgs/#Laplacians.components-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "graphAlgs",
    "title": "Laplacians.components",
    "category": "method",
    "text": "Computes the connected components of a graph. Returns them as a vector of length equal to the number of vertices. The vector numbers the components from 1 through the maximum number. For example,\n\ngr = ErdosRenyi(10,11)\nc = components(gr)\n\n10-element Array{Int64,1}:\n 1\n 1\n 1\n 1\n 2\n 1\n 1\n 1\n 3\n 2\n\n\n\n\n\n"
},

{
    "location": "graphAlgs/#Laplacians.isConnected-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "graphAlgs",
    "title": "Laplacians.isConnected",
    "category": "method",
    "text": "Returns true if graph is connected.  Calls components.\n\n\n\n\n\n"
},

{
    "location": "graphAlgs/#Laplacians.kruskal-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "graphAlgs",
    "title": "Laplacians.kruskal",
    "category": "method",
    "text": "(kruskal::SparseMatrixCSC; kind=:max) Uses Kruskal\'s algorithm to compute a minimum (or maximum) spanning tree. Set kind=:min if you want the min spanning tree. It returns it a a graph\n\n\n\n\n\n"
},

{
    "location": "graphAlgs/#Laplacians.prim-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "graphAlgs",
    "title": "Laplacians.prim",
    "category": "method",
    "text": "prim(mat::SparseMatrixCSC; kind=:max) Compute a maximum spanning tree of the matrix mat.   If kind=:min, computes a minimum spanning tree.\n\n\n\n\n\n"
},

{
    "location": "graphAlgs/#Laplacians.shortestPathTree-Tuple{Any,Any}",
    "page": "graphAlgs",
    "title": "Laplacians.shortestPathTree",
    "category": "method",
    "text": "Computes the shortest path tree, and returns it as a sparse matrix. Treats edge weights as reciprocals of lengths. For example:\n\na = [0 2 1; 2 0 3; 1 3 0]\ntr = full(shortestPathTree(sparse(a),1))\n\n3x3 Array{Float64,2}:\n 0.0  2.0  0.0\n 2.0  0.0  3.0\n 0.0  3.0  0.0\n\n\n\n\n\n"
},

{
    "location": "graphAlgs/#Laplacians.shortestPaths-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Ti}} where Ti where Tv",
    "page": "graphAlgs",
    "title": "Laplacians.shortestPaths",
    "category": "method",
    "text": "Computes the lenghts of shortest paths from start. Returns both a vector of the lenghts, and the parent array in the shortest path tree.\n\nThis algorithm treats edge weights as reciprocals of distances. DOC BETTER\n\n\n\n\n\n"
},

{
    "location": "graphAlgs/#Laplacians.vecToComps-Union{Tuple{Array{Ti,1}}, Tuple{Ti}} where Ti",
    "page": "graphAlgs",
    "title": "Laplacians.vecToComps",
    "category": "method",
    "text": "This turns a component vector, like that generated by components, into an array of arrays of indices of vertices in each component.  For example,\n\ncomps = vecToComps(c)\n\n3-element Array{Array{Int64,1},1}:\n [1,2,3,4,6,7,8]\n [5,10]\n [9]\n\n\n\n\n\n"
},

{
    "location": "graphAlgs/#Function-list-1",
    "page": "graphAlgs",
    "title": "Function list",
    "category": "section",
    "text": "Order = [:type, :function]\nPages   = [\"graphAlgs.md\"]Modules = [Laplacians]\nPages   = [\"graphAlgs.jl\"]\nPrivate = false"
},

{
    "location": "IO/#",
    "page": "IO",
    "title": "IO",
    "category": "page",
    "text": ""
},

{
    "location": "IO/#Laplacians.read_graph-Tuple{AbstractString}",
    "page": "IO",
    "title": "Laplacians.read_graph",
    "category": "method",
    "text": "adj = read_graph(fn)\n\nRead a graph from a file in IJ or IJV format. That is, each line of the file should represent one edge. Each edge should be specified by the indices of its vertices, separated by a whitespace or a comma.  If the graph is weighted, the weight should follow the second index.  For example, the unweighted complete graph on 3 vertices would appear as the file\n\n1 2\n1 3\n2 3\n\nA weighted path on 3 vertices with edge weights 1.5 and 2.5 would be\n\n1, 2, 1.5\n2, 3, 2.5\n\nThe function tries to detect the delimiter type (comma or whitespace) from the first line of the file.  The format must be consistent. Vertex indices start at 1.\n\n\n\n\n\n"
},

{
    "location": "IO/#Laplacians.writeIJV-Tuple{AbstractString,Any}",
    "page": "IO",
    "title": "Laplacians.writeIJV",
    "category": "method",
    "text": "Writes the upper portion of a matrix in ijv format, one row for each edge, separated by commas.  Only writes the upper triangular portion. The result can be read from Matlab like this:\n\n>> dl = dlmread(\'graph.txt\');\n>> a = sparse(dl(:,1),dl(:,2),dl(:,3));\n>> n = max(size(a))\n>> a(n,n) = 0;\n>> a = a + a\';\n\n\n\n\n\n"
},

{
    "location": "IO/#IO-1",
    "page": "IO",
    "title": "IO",
    "category": "section",
    "text": "Modules = [Laplacians]\nPages   = [ \"IO.jl\"]\nPrivate = false"
},

{
    "location": "solvers/#",
    "page": "solvers",
    "title": "solvers",
    "category": "page",
    "text": ""
},

{
    "location": "solvers/#Laplacians.chol_lap",
    "page": "solvers",
    "title": "Laplacians.chol_lap",
    "category": "function",
    "text": "solver = chol_lap(A::AbstractArray)\n\nUses Cholesky Factorization to solve systems in Laplacians.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.chol_sddm",
    "page": "solvers",
    "title": "Laplacians.chol_sddm",
    "category": "function",
    "text": "solveSDDM = chol_sddm(sddm::AbstractMatrix; tol, maxits, maxtime, verbose, pcgIts=Int[])\n\nThis functions wraps cholfact so that it satsfies our interface. It ignores all the keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.lapWrapSDDM-Tuple{Any,AbstractArray}",
    "page": "solvers",
    "title": "Laplacians.lapWrapSDDM",
    "category": "method",
    "text": "f = lapWrapSDDM(sddmSolver, A::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)\nf = lapWrapSDDM(sddmSolver)\n\nUses a sddmSolver to solve systems of linear equations in Laplacian matrices.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.cg",
    "page": "solvers",
    "title": "Laplacians.cg",
    "category": "function",
    "text": "x = cg(mat, b; tol, maxits, maxtime, verbose, pcgIts)\n\nsolves a symmetric linear system mat x = b.\n\nArguments\n\ntol is set to 1e-6 by default,\nmaxits defaults to Inf\nmaxtime defaults to Inf.  It measures seconds.\nverbose defaults to false\npcgIts is an array for returning the number of pcgIterations.  Default is length 0, in which case nothing is returned.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.cgLapSolver-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "solvers",
    "title": "Laplacians.cgLapSolver",
    "category": "method",
    "text": "x = cgLapSolver(A::AbstractMatrix; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])\n\nCreate a solver that uses cg to solve Laplacian systems in the laplacian of A. This just exists to satisfy our interface. It does nothing more than create the Laplacian and call cg on each connected component.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.cgSolver",
    "page": "solvers",
    "title": "Laplacians.cgSolver",
    "category": "function",
    "text": "x = cgSolver(mat; tol, maxits, maxtime, verbose, pcgIts)\n\ncreates a solver for a PSD system mat. The parameters are as described in cg.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.pcg",
    "page": "solvers",
    "title": "Laplacians.pcg",
    "category": "function",
    "text": "x = pcg(mat, b, pre; tol, maxits, maxtime, verbose, pcgIts, stag_test)`\n\nsolves a symmetric linear system using preconditioner pre.\n\nArguments\n\npre can be a function or a matrix.  If a matrix, a function to solve it is created with cholFact.\ntol is set to 1e-6 by default,\nmaxits defaults to Inf\nmaxtime defaults to Inf.  It measures seconds.\nverbose defaults to false\npcgIts is an array for returning the number of pcgIterations.  Default is length 0, in which case nothing is returned.\nstag_test=k stops the code if rho[it] > (1-1/k) rho[it-k].  Set to 0 to deactivate.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.pcgLapSolver-Tuple{AbstractArray{T,2} where T,AbstractArray{T,2} where T}",
    "page": "solvers",
    "title": "Laplacians.pcgLapSolver",
    "category": "method",
    "text": "x = pcgLapSolver(A, B; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[])\n\nCreate a solver that uses pcg to solve Laplacian systems in A Specialized for the case when the preconditioner the Laplacian matrix of B. It solves the preconditioner by Cholesky Factorization.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.pcgSolver",
    "page": "solvers",
    "title": "Laplacians.pcgSolver",
    "category": "function",
    "text": "x = pcgSolver(mat, pre; tol, maxits, maxtime, verbose, pcgIts)\n\ncreates a solver for a PSD system using preconditioner pre. The parameters are as described in pcg.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.ApproxCholParams",
    "page": "solvers",
    "title": "Laplacians.ApproxCholParams",
    "category": "type",
    "text": "params = ApproxCholParams(order, output)\n\norder can be one of\n\n:deg (by degree, adaptive),\n:wdeg (by original wted degree, nonadaptive),\n:given\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.approxchol_lap-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "solvers",
    "title": "Laplacians.approxchol_lap",
    "category": "method",
    "text": "solver = approxchol_lap(a); x = solver(b);\nsolver = approxchol_lap(a; tol::Real=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())\n\nA heuristic by Daniel Spielman inspired by the linear system solver in https://arxiv.org/abs/1605.02353 by Rasmus Kyng and Sushant Sachdeva.  Whereas that paper eliminates vertices one at a time, this eliminates edges one at a time.  It is probably possible to analyze it. The ApproxCholParams let you choose one of three orderings to perform the elimination.\n\nApproxCholParams(:given) - in the order given.   This is the fastest for construction the preconditioner, but the slowest solve.\nApproxCholParams(:deg) - always eliminate the node of lowest degree.   This is the slowest build, but the fastest solve.\nApproxCholParams(:wdeg) - go by a perturbed order of wted degree.\n\nFor more info, see http://danspielman.github.io/Laplacians.jl/latest/usingSolvers/index.html\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.approxchol_sddm",
    "page": "solvers",
    "title": "Laplacians.approxchol_sddm",
    "category": "function",
    "text": "solver = approxchol_sddm(sddm); x = solver(b);\nsolver = approxchol_sddm(sddm; tol=1e-6, maxits=1000, maxtime=Inf, verbose=false, pcgIts=Int[], params=ApproxCholParams())\n\nSolves sddm systems by wrapping approxchol_lap. Not yet optimized directly for sddm.\n\nFor more info, see http://danspielman.github.io/Laplacians.jl/latest/usingSolvers/index.html\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.condNumber-Tuple{Any,Any}",
    "page": "solvers",
    "title": "Laplacians.condNumber",
    "category": "method",
    "text": "cn = condNumber(a, ldli; verbose=false)\n\nGiven an adjacency matrix a and an ldli computed by approxChol, this computes the condition number.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.augTreeLap-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "solvers",
    "title": "Laplacians.augTreeLap",
    "category": "method",
    "text": "solver = augTreeLap(A; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params=AugTreeParams())\n\nAn \"augmented spanning tree\" solver for Laplacians.  It works by adding edges to a low stretch spanning tree.  It calls augTreePrecon to form the preconditioner.  params has entries\n\nparams.treeAlg default to akpw\nparams.opt if true, it interacts with cholmod to choose a good number of edges to add back.  If false, it adds back 2*sqrt(n).\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.augTreeLapPrecon-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "solvers",
    "title": "Laplacians.augTreeLapPrecon",
    "category": "method",
    "text": "pre = augTreeLapPrecon{Tv,Ti}(A; params=AugTreeParams())\n\nThis is an augmented spanning tree preconditioner for Laplacians. It takes as optional input a tree growing algorithm. It adds back 2sqrt(n) edges via augmentTree: the sqrt(n) of highest stretch and another sqrt(n) sampled according to stretch. For most purposes, one should directly call augTreeLapSolver.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.augTreePrecon-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "solvers",
    "title": "Laplacians.augTreePrecon",
    "category": "method",
    "text": "pre = augTreePrecon{Tv,Ti}(ddmat::SparseMatrixCSC{Tv,Ti}; params=AugTreeParams())\n\nThis is an augmented spanning tree preconditioner for diagonally dominant linear systems.  It takes as optional input a tree growing algorithm. It adds back 2sqrt(n) edges via augmentTree: the sqrt(n) of highest stretch and another sqrt(n) sampled according to stretch. For most purposes, one should directly call augTreeSolver.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.augTreeSddm-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "solvers",
    "title": "Laplacians.augTreeSddm",
    "category": "method",
    "text": "solver = augTreeSddm(sddm; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[],  params=AugTreeParams())\n\nAn \"augmented spanning tree\" solver for positive definite diagonally dominant matrices.  It works by adding edges to a low stretch spanning tree.  It calls augTreePrecon to form the preconditioner.  params has entries\n\nparams.treeAlg default to akpw\nparams.opt if true, it interacts with cholmod to choose a good number of edges to add back.  If false, it adds back 2*sqrt(n).\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.augmentTree-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},SparseArrays.SparseMatrixCSC{Tv,Ti},Ti}} where Ti where Tv",
    "page": "solvers",
    "title": "Laplacians.augmentTree",
    "category": "method",
    "text": "B = augmentTree{Tv,Ti}(tree, A, k)\n\nTakes as input a tree and an adjacency matrix of a graph. It then computes the stretch of every edge of the graph wrt the tree.  It then adds back the k edges of highest stretch, and k edges sampled according to stretch.\n\nThis is the old alg.  We now recommend using augmentTreeOpt.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.KMPParams",
    "page": "solvers",
    "title": "Laplacians.KMPParams",
    "category": "type",
    "text": "Parameters for the KMP solver\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.KMPLapSolver-Tuple{Any}",
    "page": "solvers",
    "title": "Laplacians.KMPLapSolver",
    "category": "method",
    "text": "lapSolver = KMPLapSolver(A; verbose, tol, maxits, maxtime, pcgIts, params::KMPParams)\n\nSolves linear equations in the Laplacian of graph with adjacency matrix A.\n\nBased on the paper \"Approaching optimality for solving SDD systems\" by Koutis, Miller, and Peng, <i>SIAM Journal on Computing</i>, 2014.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.KMPSDDMSolver-Tuple{Any}",
    "page": "solvers",
    "title": "Laplacians.KMPSDDMSolver",
    "category": "method",
    "text": "sddmSolver = KMPSDDMSolver(mat; verbose, tol, maxits, maxtime, pcgIts, params::KMPParams)\n\nSolves linear equations in symmetric, diagonally dominant matrices with non-positive off-diagonals.  Based on the paper \"Approaching optimality for solving SDD systems\" by Koutis, Miller, and Peng, <i>SIAM Journal on Computing</i>, 2014.\n\n\n\n\n\n"
},

{
    "location": "solvers/#Laplacians.harmonic_interp-Tuple{Any,Array{T,1} where T,Array{T,1} where T}",
    "page": "solvers",
    "title": "Laplacians.harmonic_interp",
    "category": "method",
    "text": "x = harmonic_interp(adj_mat, S, vals; tol=1e-6)\n\nInterpolates a function on a graph, given by its adjacency matrix, by minizing the Laplacian quadratic form subject to the boundary conditions that x[S[i]] = vals[i] for i in S.\n\nThis is the algorithm sometimes known as Label Propagation, or Semi-Supervised Learning on Graphs.  The idea comes from the paper \"Semi-Supervised Learning Using Gaussian Fields and Harmonic Functions\" by Zhu, Gharamani, and Lafferty from ICML 2003.\n\nThis version might fail for disconnected graphs. You can check if a graph is connected with isConnected(adj_mat).\n\n\n\n\n\n"
},

{
    "location": "solvers/#Linear-Equation-Solvers-1",
    "page": "solvers",
    "title": "Linear Equation Solvers",
    "category": "section",
    "text": "For more, see the page on using solvers.Order = [:type, :function]\nPages   = [\"solvers.md\"]Modules = [Laplacians]\nPages   = [\"solverInterface.jl\", \"pcg.jl\",\"approxChol.jl\", \"augTreeSolver.jl\", \"samplingSolver.jl\",\"KMPSolver.jl\", \"externalSolvers.jl\", \"harmonic.jl\"]\nPrivate = false"
},

{
    "location": "sparsification/#",
    "page": "sparsification",
    "title": "sparsification",
    "category": "page",
    "text": ""
},

{
    "location": "sparsification/#Laplacians.sparsify-Tuple{Any}",
    "page": "sparsification",
    "title": "Laplacians.sparsify",
    "category": "method",
    "text": "as = sparsify(a; ep=0.5)\n\nApply Spielman-Srivastava sparsification: sampling by effective resistances. ep should be less than 1.\n\n\n\n\n\n"
},

{
    "location": "sparsification/#Laplacians.approxQual-Tuple{Any,Any}",
    "page": "sparsification",
    "title": "Laplacians.approxQual",
    "category": "method",
    "text": "eps = approxQual(graph1, graph2; tol=1e-5, verbose=false)\n\nComputes the eps for which graph1 and graph2 are eps approximations of each other. That is, L1 <= (1+eps) L2, and vice versa.\n\nIt is randomized, so you might want to run it again if you don\'t trust the answers.\n\n\n\n\n\n"
},

{
    "location": "sparsification/#Laplacians.conditionNumber-Tuple{SparseArrays.SparseMatrixCSC,Function}",
    "page": "sparsification",
    "title": "Laplacians.conditionNumber",
    "category": "method",
    "text": "kappa = conditionNumber(graph, precon; tol=1e-5, verbose=false)\n\nComputes the relative condition number of graph and a preconditioning function.\n\nIt is randomized, so you might want to run it again if you don\'t trust the answers.\n\n\n\n\n\n"
},

{
    "location": "sparsification/#Laplacians.conditionNumber-Tuple{SparseArrays.SparseMatrixCSC,SparseArrays.SparseMatrixCSC}",
    "page": "sparsification",
    "title": "Laplacians.conditionNumber",
    "category": "method",
    "text": "kapps = conditionNumber(graph1, graph2; tol=1e-5, verbose=false)\n\nComputes the relative condition number of graph1 and graph2.\n\nIt is randomized, so you might want to run it again if you don\'t trust the answers.\n\n\n\n\n\n"
},

{
    "location": "sparsification/#Laplacians.support-Tuple{Any,Any}",
    "page": "sparsification",
    "title": "Laplacians.support",
    "category": "method",
    "text": "sup12, sup21 = support(graph1, graph2; tol=1e-5)\n\nComputes the support of graph1 wrt graph2, and the other way around. It is randomized, so you might want to run it again if you don\'t trust the answers.\n\n\n\n\n\n"
},

{
    "location": "sparsification/#sparsification-1",
    "page": "sparsification",
    "title": "sparsification",
    "category": "section",
    "text": "Order = [:type, :function]\nPages   = [\"sparsification.md\"]Modules = [Laplacians]\nPages   = [\"sparsify.jl\",\"conditionNumber.jl\"]\nPrivate = false"
},

{
    "location": "akpw/#",
    "page": "akpw",
    "title": "akpw",
    "category": "page",
    "text": ""
},

{
    "location": "akpw/#Laplacians.akpw-Tuple{Any}",
    "page": "akpw",
    "title": "Laplacians.akpw",
    "category": "method",
    "text": "tree = akpw(graph; ver=0)\n\nComputes a low stretch spanning tree of graph, and returns it as a graph. The default version is 0.  In event of emergency, one can try ver=2.  It is usually slower, but might have slightly better stretch.\n\n\n\n\n\n"
},

{
    "location": "akpw/#Laplacians.akpwU-Tuple{Any}",
    "page": "akpw",
    "title": "Laplacians.akpwU",
    "category": "method",
    "text": "tree = akpwU(graph)\n\nComputes a low stretch spanning tree of an unweighted graph, and returns it as a graph.\n\n\n\n\n\n"
},

{
    "location": "akpw/#AKPW-1",
    "page": "akpw",
    "title": "AKPW",
    "category": "section",
    "text": "Also see the page on Low Stretch Spanning TreesOrder = [:type, :function]\nPages   = [\"akpw.md\"]Modules = [Laplacians]\nPages   = [\"akpw.jl\"]\nPrivate = false"
},

{
    "location": "treeAlgs/#",
    "page": "treeAlgs",
    "title": "treeAlgs",
    "category": "page",
    "text": ""
},

{
    "location": "treeAlgs/#Laplacians.comp_stretches-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},SparseArrays.SparseMatrixCSC{Tv,Ti}}} where Ti where Tv",
    "page": "treeAlgs",
    "title": "Laplacians.comp_stretches",
    "category": "method",
    "text": "Compute the stretched of every edge in mat with respect to the tree tree. Returns the answer as a sparse matrix with the same nonzero structure as mat. Assumes that mat is symmetric. tree should be the adjacency matrix of a spanning tree.\n\n\n\n\n\n"
},

{
    "location": "treeAlgs/#Tree-Algorithms-1",
    "page": "treeAlgs",
    "title": "Tree Algorithms",
    "category": "section",
    "text": "Order = [:type, :function]\nPages   = [\"treeAlgs.md\"]Modules = [Laplacians]\nPages   = [\"treeAlgs.jl\"]\nPrivate = false"
},

{
    "location": "randTrees/#",
    "page": "randTrees",
    "title": "randTrees",
    "category": "page",
    "text": ""
},

{
    "location": "randTrees/#Laplacians.randishKruskal-Tuple{SparseArrays.SparseMatrixCSC}",
    "page": "randTrees",
    "title": "Laplacians.randishKruskal",
    "category": "method",
    "text": "tree = randishKruskal(A)\n\nA heuristic for computing low-stretch spanning trees.  Where Kruskal\'s MST algorithm adds edges in order of weight, this algorithm adds them at random with probability proportional to their weight.\n\n\n\n\n\n"
},

{
    "location": "randTrees/#Laplacians.randishPrim-Union{Tuple{SparseArrays.SparseMatrixCSC{Tval,Tind}}, Tuple{Tind}, Tuple{Tval}} where Tind where Tval",
    "page": "randTrees",
    "title": "Laplacians.randishPrim",
    "category": "method",
    "text": "tree = randishPrim(A)\n\nA heuristic for computing low-stretch spanning trees.  Where Prim\'s MST algorithm grows a cluster by always adding the edge on the boundary of maximum weight, this algorithm adds a boundary edge with probability proportional to its weight.\n\n\n\n\n\n"
},

{
    "location": "randTrees/#randTrees-1",
    "page": "randTrees",
    "title": "randTrees",
    "category": "section",
    "text": "Order = [:type, :function]\nPages   = [\"randTrees.md\"]Modules = [Laplacians]\nPages   = [\"randTrees.jl\"]\nPrivate = false"
},

{
    "location": "localClustering/#",
    "page": "localClustering",
    "title": "localClustering",
    "category": "page",
    "text": ""
},

{
    "location": "localClustering/#Laplacians.dumbRefineCut-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1}}} where Ti where Tv",
    "page": "localClustering",
    "title": "Laplacians.dumbRefineCut",
    "category": "method",
    "text": "Modify a cluster by passing through all the vertices exactly once and \nadding/removing them based on the value of (Deg_external - Deg_Internal).\n\n\n\n\n\n"
},

{
    "location": "localClustering/#Laplacians.refineCut-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1}}} where Ti where Tv",
    "page": "localClustering",
    "title": "Laplacians.refineCut",
    "category": "method",
    "text": "Modify a cluster by adding or removing vertices by picking at each step \nthe vertex that has the maximum value of (Deg_external - Deg_Internal).\nEach vertex can be added in/removed only once.\n\n\n\n\n\n"
},

{
    "location": "localClustering/#Laplacians.apr-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1},Float64,Float64}} where Ti where Tv",
    "page": "localClustering",
    "title": "Laplacians.apr",
    "category": "method",
    "text": "Computes an approximate page rank vector from a starting set s, an alpha and an epsilon The algorithm follows the Anderson,Chung,Lang paper and Dan Spielman\'s lecture notes\n\n\n\n\n\n"
},

{
    "location": "localClustering/#Laplacians.localImprove-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1}}} where Ti where Tv",
    "page": "localClustering",
    "title": "Laplacians.localImprove",
    "category": "method",
    "text": "localImprove{Tv,Ti}(G::SparseMatrixCSC{Tv,Ti}, A::Array{Int64,1}; epsSigma=-1.0, err=1e-10, maxSize = max(G.n, G.m)\n\nThe LocalImprove function, from the Orrechia-Zhu paper. Given a graph and an initial set, finds a set of smaller conductance based on the starting set using a localized version of max-flow.\n\nSmall discussion: When adding in the neighbors of the initial component, if the resulting  conductance is worse than the initial one,  the algorithm will add more and more vertices until hitting a better conductance. However, if we fix a certain  maximum size for our component,  it might be the case that this new conductance will always be worse than what we had initially. Thus, if we run the algorithm with a small maxSize,  our initial conductance might be the best solution we can raech.\n\nG is the given graph, A is the initial set \nepsSigma is a measure of the quality of the returning set (the smaller the better). It\'s defaulted to volume(A) / volume(V - A)\nerr is the numerical error considered throughout the algorithm. It\'s defaulted to 1e-10\nmaxSize is the maximum allowed size for the flow graph at any iteration of the algorithm. It\'s defaulted to |V|\n\n\n\n\n\n"
},

{
    "location": "localClustering/#Laplacians.prn-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1},Float64,Int64}} where Ti where Tv",
    "page": "localClustering",
    "title": "Laplacians.prn",
    "category": "method",
    "text": "prn{Tv, Ti}(G::SparseMatrixCSC{Tv,Ti}, s::Array{Int64,1}, phi::Float64, b::Int64)\n\nThe PageRank-Nibble cutting algorithm from the Anderson/Chung/Lang paper\n\ns is a set of starting vertices, phi is a constant in (0, 1], and b is an integer in [1, [log m]]\n\nphi is a bound on the quality of the conductance of the cut - the smaller the phi, the higher the quality.  b is used to handle precision throughout the algorithm - the higher the b, the greater the precision.\n\n\n\n\n\n"
},

{
    "location": "localClustering/#Local-Clustering-1",
    "page": "localClustering",
    "title": "Local Clustering",
    "category": "section",
    "text": "This is a collection of clustering related algorithms,   based on Approximate Personal PageRank, and improvement by local   flow computations.   It needs more documentation.For now, see the Local Clustering NotebookPages   = [\"localClustering.md\"]Modules = [Laplacians]\nPages   = [\"cutHeuristics.jl\",\"localClustering.jl\"]"
},

{
    "location": "privateFuncs/#",
    "page": "Private Functions",
    "title": "Private Functions",
    "category": "page",
    "text": ""
},

{
    "location": "privateFuncs/#Laplacians.ApproxCholPQ",
    "page": "Private Functions",
    "title": "Laplacians.ApproxCholPQ",
    "category": "type",
    "text": "An approximate priority queue.   Items are bundled together into doubly-linked lists with all approximately the same key.   minlist is the min list we know to be non-empty.   It should always be a lower bound.   keyMap maps keys to lists\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.IJV-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.IJV",
    "category": "method",
    "text": "ijv = IJV(A::SparseMatrixCSC)\n\nConvert a sparse matrix to an IJV.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.LDLinv",
    "page": "Private Functions",
    "title": "Laplacians.LDLinv",
    "category": "type",
    "text": "LDLinv contains the information needed to solve the Laplacian systems.   It does it by applying Linv, then Dinv, then Linv (transpose).   But, it is specially constructed for this particular solver.   It does not explicitly make the matrix triangular.   Rather, col[i] is the name of the ith col to be eliminated\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.LLmatp",
    "page": "Private Functions",
    "title": "Laplacians.LLmatp",
    "category": "type",
    "text": "LLmatp is the data structure used to maintain the matrix during elimination.   It stores the elements in each column in a singly linked list (only next ptrs)   Each element is an LLp (linked list pointer).   The head of each column is pointed to by cols.\n\nWe probably can get rid of degs - as it is only used to store initial degrees.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.LLp",
    "page": "Private Functions",
    "title": "Laplacians.LLp",
    "category": "type",
    "text": "LLp elements are all in the same column.   row tells us the row, and val is the entry.   val is set to zero for some edges that we should remove.   next gives the next in the column.  It points to itself to terminate.   reverse is the index into lles of the other copy of this edge,   since every edge is stored twice as we do not know the order of elimination in advance.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.addToGPrime-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Array{Tuple{Int64,Float64},1},1},Dict{Int64,Int64},Dict{Int64,Int64},Int64,Int64,Float64,Float64,Int64}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.addToGPrime",
    "category": "method",
    "text": "Add a new vertex to GPrime \n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.approxCholPQDec!-Union{Tuple{Tind}, Tuple{Laplacians.ApproxCholPQ{Tind},Any}} where Tind",
    "page": "Private Functions",
    "title": "Laplacians.approxCholPQDec!",
    "category": "method",
    "text": "Decrement the key of element i\nThis could crash if i exceeds the maxkey\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.approxCholPQInc!-Union{Tuple{Tind}, Tuple{Laplacians.ApproxCholPQ{Tind},Any}} where Tind",
    "page": "Private Functions",
    "title": "Laplacians.approxCholPQInc!",
    "category": "method",
    "text": "Increment the key of element i\nThis could crash if i exceeds the maxkey\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.approxchol_lapChol-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.approxchol_lapChol",
    "category": "method",
    "text": "This variation of approxChol creates a cholesky factor to do the elimination. It has not yet been optimized, and does not yet make the cholesky factor lower triangular\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.augmentTreeOpt-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},SparseArrays.SparseMatrixCSC{Tv,Ti}}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.augmentTreeOpt",
    "category": "method",
    "text": "B = augmentTreeOpt{Tv,Ti}(tree, A, params)\n\nTakes as input a tree and an adjacency matrix of a graph. It then computes the stretch of every edge of the graph wrt the tree.  It uses cholmod to decide how many edge to add back, shooting for nnzLfac times n entries in the factored augmented tree, with a number of flops to factor equal to nnz(a)*flopsfac. The edges to add back are then choen at random.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.blockSolver-Tuple{Any,Any}",
    "page": "Private Functions",
    "title": "Laplacians.blockSolver",
    "category": "method",
    "text": "Apply the ith solver on the ith component\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.complete_bipartite_graph_ijv-Tuple{Integer}",
    "page": "Private Functions",
    "title": "Laplacians.complete_bipartite_graph_ijv",
    "category": "method",
    "text": "ijv = complete_bipartite_graph_ijv(n)\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.complete_graph_ijv-Tuple{Integer}",
    "page": "Private Functions",
    "title": "Laplacians.complete_graph_ijv",
    "category": "method",
    "text": "ijv = complete_graph_ijv(n)\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.empty_graph_ijv-Tuple{Integer}",
    "page": "Private Functions",
    "title": "Laplacians.empty_graph_ijv",
    "category": "method",
    "text": "ijv = empty_graph_ijv(n)\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.extendMatrix-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Tv,1}}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.extendMatrix",
    "category": "method",
    "text": "Add a new vertex to a with weights to the other vertices corresponding to diagonal surplus weight.\n\nThis is an efficient way of writing [a d; d\' 0]\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.firstn-Tuple{Laplacians.IJV,Integer}",
    "page": "Private Functions",
    "title": "Laplacians.firstn",
    "category": "method",
    "text": "b = firstn(a::IJV, n::Integer)\n\nOnly keep the first n vertices of a.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.forceLap-Tuple{AbstractArray}",
    "page": "Private Functions",
    "title": "Laplacians.forceLap",
    "category": "method",
    "text": "la = forceLap(a)\n\nCreate a Laplacian matrix from an adjacency matrix. If the input looks like a Laplacian, throw a warning and convert it.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.generalizedNecklace-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},SparseArrays.SparseMatrixCSC,Int64}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.generalizedNecklace",
    "category": "method",
    "text": "graph = generalizedNecklace(A, H, k::Int64)\n\nConstructs a generalized necklace graph starting with two graphs A and H. The resulting new graph will be constructed by expanding each vertex in H to an instance of A. k random edges will be generated between components. Thus, the resulting graph may have weighted edges.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.getCutSet-Tuple{Array{Array{Tuple{Int64,Float64},1},1},Int64,Int64}",
    "page": "Private Functions",
    "title": "Laplacians.getCutSet",
    "category": "method",
    "text": "Get the min cut from the source - return all vertices in the cut besides the source \n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.initDictCol!-Tuple{Any,Any,Any}",
    "page": "Private Functions",
    "title": "Laplacians.initDictCol!",
    "category": "method",
    "text": "initDictCol!(dic, name, typ)\n\nFor a dictionary in which each key indexes an array. If dic does not contain an entry of name, create with set to Array(typ,0).\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.initGPrime-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1},Dict{Int64,Int64},Dict{Int64,Int64},Float64,Int64}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.initGPrime",
    "category": "method",
    "text": "Initialize GPrime with the set A and edges of type s->u\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.lapWrapComponents-Tuple{Any,AbstractArray}",
    "page": "Private Functions",
    "title": "Laplacians.lapWrapComponents",
    "category": "method",
    "text": "f = lapWrapComponents(solver, a::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)\n\nApplies a Laplacian solver that satisfies our interface to each connected component of the graph with adjacency matrix a. Passes kwargs on the solver.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.lapWrapConnected-Tuple{Any,AbstractArray{T,2} where T}",
    "page": "Private Functions",
    "title": "Laplacians.lapWrapConnected",
    "category": "method",
    "text": "f = lapWrapConnected(sddmSolver, a::AbstractMatrix; kwargs...)\n\nApplies a sddmSolver to the Laplacian of the adjacency matrix a of a connected graph. Passes on kwargs to the solver. sddmSolver should be a solver that obeys the interface.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.ldli2Chol-Tuple{Any}",
    "page": "Private Functions",
    "title": "Laplacians.ldli2Chol",
    "category": "method",
    "text": "L = ldli2Chol(ldli)\n\nThis produces a matrix L so that L L^T approximate the original Laplacians. It is not quite a Cholesky factor, because it is off by a perm (and the all-1s vector orthogonality.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.localBlockFlow-Tuple{Array{Array{Tuple{Int64,Float64},1},1},Int64,Int64}",
    "page": "Private Functions",
    "title": "Laplacians.localBlockFlow",
    "category": "method",
    "text": "Compute block flow between s and t\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.localFlow-Union{Tuple{Ti}, Tuple{Tv}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1},Float64,Float64}, Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti},Array{Int64,1},Float64,Float64,Any}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.localFlow",
    "category": "method",
    "text": "The LocalFlow function, from the Orecchia-Zhu paper \n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.path_graph_ijv-Tuple{Integer}",
    "page": "Private Functions",
    "title": "Laplacians.path_graph_ijv",
    "category": "method",
    "text": "ijv = path_graph_ijv(n::Int64)\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.print_ll_col-Tuple{Laplacians.LLMatOrd,Int64}",
    "page": "Private Functions",
    "title": "Laplacians.print_ll_col",
    "category": "method",
    "text": "Print a column in an LLMatOrd matrix.   This is here for diagnostics.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.print_ll_col-Tuple{Laplacians.LLmatp,Int64}",
    "page": "Private Functions",
    "title": "Laplacians.print_ll_col",
    "category": "method",
    "text": "Print a column in an LLmatp matrix.   This is here for diagnostics.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.pure_random_ijv-Tuple{Integer}",
    "page": "Private Functions",
    "title": "Laplacians.pure_random_ijv",
    "category": "method",
    "text": "a = pure_random_ijv(n::Integer; verbose=false, prefix=\"\", ver=Vcur)\n\nChooses among pathgraph, ringgraph, gridgraph, completebinarytree, randgenring, growngraph and ErdosRenyiClusterFix. It can produce a disconnected graph. For code that always produces a connected graph (and is the same as with Julia v0.6, use purerandomijv_v6)\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.pushSpeedResult!-Tuple{Any,Any,Any}",
    "page": "Private Functions",
    "title": "Laplacians.pushSpeedResult!",
    "category": "method",
    "text": "ret is the answer returned by a speed test. This pushed it into the dictionary on which we are storing the tests.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.rand_regular_bipartite-Tuple{Any,Any}",
    "page": "Private Functions",
    "title": "Laplacians.rand_regular_bipartite",
    "category": "method",
    "text": "a = rand_regular_bipartite(n,k)\n\nRandom k-regular bipartite graph between two sets of n vertices. No repeat edges, so can take a long time to build of k is close to n.\n\nReturns a (possibly) asymmetric matrix that contains the upper-right block.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.ring_graph_ijv-Tuple{Integer}",
    "page": "Private Functions",
    "title": "Laplacians.ring_graph_ijv",
    "category": "method",
    "text": "ijv = ring_graph_ijv(n)\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.sampleByWeight-Tuple{Any}",
    "page": "Private Functions",
    "title": "Laplacians.sampleByWeight",
    "category": "method",
    "text": "ind = sampleByWeight(wt; ver=Vcur)\n\nsample an index with probability proportional to its weight given here\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.sddmWrapLap-Tuple{Any,AbstractArray}",
    "page": "Private Functions",
    "title": "Laplacians.sddmWrapLap",
    "category": "method",
    "text": "f = sddmWrapLap(lapSolver, sddm::AbstractArray; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)\n\nUses a lapSolver to solve systems of linear equations in sddm matrices.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.sortSet-Union{Tuple{Ti}, Tuple{Array{Ti,1},Ti}} where Ti",
    "page": "Private Functions",
    "title": "Laplacians.sortSet",
    "category": "method",
    "text": "Given a set of integers, set between 1 and n, return a sorted version of them\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.star_graph_ijv-Tuple{Integer}",
    "page": "Private Functions",
    "title": "Laplacians.star_graph_ijv",
    "category": "method",
    "text": "ijv = star_graph_ijv(n::Int64)\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.testZeroDiag-Tuple{Any}",
    "page": "Private Functions",
    "title": "Laplacians.testZeroDiag",
    "category": "method",
    "text": "testZeroDiag(a)\n\nReturns true if a has zero diagonal, false otherwise\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.treeDepthDFS-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv,Ti}}, Tuple{Ti}, Tuple{Tv}} where Ti where Tv",
    "page": "Private Functions",
    "title": "Laplacians.treeDepthDFS",
    "category": "method",
    "text": "Compute the vector of depths in a tree that is in DFS order, with the root at the first position, and the leaves at the end\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.uniformWeight_ver-Tuple{Type{Val{6}},SparseArrays.SparseMatrixCSC}",
    "page": "Private Functions",
    "title": "Laplacians.uniformWeight_ver",
    "category": "method",
    "text": "wted = uniformWeight(unwted)\n\nPut a uniform [0,1] weight on every edge.  This is an example of how to use mapweight.\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.wrapCapture-Tuple{Function,Any,Any}",
    "page": "Private Functions",
    "title": "Laplacians.wrapCapture",
    "category": "method",
    "text": "f = wrapCapture(solver::Function, mats, rhss; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)\n\nThis wraps a solver so that we can capture all the matrices that it solves and all the right-hand-sides. Those are pushed into the arrays mats and rhss. For example\n\njulia> mats = []\njulia> rhss = []\njulia> solver = wrapCapture(approxchol_lap, mats, rhss)\njulia> a = chimera(10)\njulia> f = solver(a);\njulia> size(mats[1])\n(10,10)\njulia> b = randn(10)\njulia> x = f(b);\njulia> rhss\n1-element Array{Any,1}:\n [0.404962,-0.827718,0.704616,-0.403223,0.204891,-0.505589,0.907015,1.90266,-0.438115,0.0464351]\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.wrapCaptureRhs-Tuple{Function,Any}",
    "page": "Private Functions",
    "title": "Laplacians.wrapCaptureRhs",
    "category": "method",
    "text": "f = wrapCaptureRhs(sola::Function, rhss; tol::Real=1e-6, maxits=Inf, maxtime=Inf, verbose=false, pcgIts=Int[], params...)\n\nCaptures all the right-hand-sides that are passed to the solver sola.  It pushes them into an array called rhhs. For example\n\njulia> rhss = []\njulia> a = wted_chimera(100)\njulia> sola = approxchol_lap(a)\njulia> wrappedSolver = wrapCaptureRhs(sola,rhss)\njulia> b = randn(100)\njulia> x = wrappedSolver(b,verbose=true)\n\nPCG BLAS stopped after: 0.0 seconds and 11 iterations with relative error 3.160275810360986e-7.\n\njulia> length(rhss[1])\n\n100\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Laplacians.wrapInterface-Tuple{Function,AbstractArray{T,2} where T}",
    "page": "Private Functions",
    "title": "Laplacians.wrapInterface",
    "category": "method",
    "text": "solveA = wrapInterface(solver::Function, A::AbstractMatrix; tol, maxits, maxtime, verbose, pcgIts=Int[],params...)\nsolverConstructor = wrapInterface(A::AbstractMatrix; tol, maxits, maxtime, verbose, pcgIts=Int[],params...)\n\nReturns a function that discards tol, maxits, maxtime and verbose, sets pcgIts to 0 (because it might not be using pcg), and passes whatever params are left to the solver.\n\nExamples\n\njulia> a = randn(5,5);\njulia> a = a * a\';\njulia> solvea = wrapInterface(X->cholesky(X,Val(true)), a, maxits=100, verbose=true);\njulia> b = randn(5,1);\njulia> norm(a*solvea(b, verbose=false)-b)\n1.575705319704736e-14\n\njulia> f = wrapInterface(X->cholesky(X,Val(true)))\njulia> solvea = f(a, maxits=1000, maxtime = 1)\njulia> norm(a*solvea(b, verbose=false, maxtime = 10)-b)\n1.575705319704736e-14\n\n\n\n\n\n"
},

{
    "location": "privateFuncs/#Unexported-(Private)-functions.-1",
    "page": "Private Functions",
    "title": "Unexported (Private) functions.",
    "category": "section",
    "text": "This is a list of all unexported functions and types from Laplacians.Pages   = [\"privateFuncs.md\"]Public = false\nModules = [Laplacians]"
},

{
    "location": "indexOfAll/#",
    "page": "All of the above",
    "title": "All of the above",
    "category": "page",
    "text": ""
},

{
    "location": "indexOfAll/#Index-of-all-exported-1",
    "page": "All of the above",
    "title": "Index of all exported",
    "category": "section",
    "text": "This is an index of all the exported methods. We would include the docstrings, but Documenter.jl does not let us.modules = [Laplacians]"
},

]}
