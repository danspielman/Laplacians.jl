"""
    e = fiedler(A;  nev = 1, tol = 0.0 ) 

Compute the Fiedler vector and value of the graph with adjacency matrix `A`.
Computes them using `eigs` from Arpack.jl, and returns `e` in the same format.
So, `e[1]` is the fiedler value, and `e[2]` is the vector.
If `nev` is more than 1, then it returns `nev` values and vectors.

Note that for small matrices, or for well-conditioned ones, it will be faster to simply use
~~~    
    e = eigs(lap(a), which = :SM, nev = 1)
~~~
"""
function fiedler(a::AbstractArray; nev = 1, tol = 0.0)
    f = approxchol_lap(a) 
    op = Laplacians.SqLinOp(true,1.0,size(a,1),f)
    e = eigs(op, which=:LM, nev=nev, tol=tol)
    e[1] .= 1 ./ e[1]
    return e
end
