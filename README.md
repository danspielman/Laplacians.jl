# Laplacians.jl
Laplacians is a package containing graph algorithms, with an emphsasis on tasks related to spectral and algebraic graph theory. It contains (and will contain more) code for solving systems of linear equations in graph Laplacians, low stretch spanning trees, sparsifiation, clustering, local clustering, and optimization on graphs.

All graphs are represented by sparse adjacency matrices. This is both for speed, and because our main concerns are algebraic tasks. It does not handle dynamic graphs. It would be very slow to implement dynamic graphs this way.

The documentation may be found in
[site/about/index.html](site/about/index.html).
If you are using a Mac, you can access the documentation by typing

~~~
open site/about/index.html
~~~

