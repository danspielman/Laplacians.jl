#=

Code for representing sparse matrices in IJV format.
It is often more efficient to manipulate these.
Our graph generators with _ijv at the end of their name produce these.

=#

import Base: ==, hash, +, *, transpose, adjoint

mutable struct IJV{Tv,Ti}
    n::Ti
    nnz::Ti
    i::Array{Ti,1}
    j::Array{Ti,1}
    v::Array{Tv,1}
end

==(a::IJV, b::IJV) = 
    a.n == b.n &&
    a.nnz == b.nnz &&
    a.i == b.i &&
    a.j == b.j &&
    a.v == b.v

hash(a::IJV) = 
    hash((a.n, a.nnz, a.i, a.j, a.v), hash(IJV))

hash(a::IJV, h::UInt) = hash(hash(a), h)

function +(a::IJV, b::IJV)
    @assert a.n == b.n
    IJV(a.n, a.nnz+b.nnz,
        [a.i; b.i],
        [a.j; b.j],
        [a.v; b.v])
end

function *(a::IJV, x::Number)
    ijv = deepcopy(a)
    ijv.v .*= x

    return ijv
end

*(x::Number, a::IJV) = *(a, x)

transpose(ijv::IJV) = IJV(ijv.n, ijv.nnz, ijv.j, ijv.i, ijv.v)
adjoint(ijv::IJV) = IJV(ijv.n, ijv.nnz, ijv.j, ijv.i, adjoint.(ijv.v))


"""
    ijv = IJV(A::SparseMatrixCSC)
Convert a sparse matrix to an IJV.
"""
function IJV(A::SparseMatrixCSC{Tv, Ti}) where {Tv, Ti}
    (ai,aj,av) = findnz(A)
    IJV{Tv,Ti}(A.n, nnz(A), ai, aj, av)
end


import SparseArrays.sparse

sparse(ijv::IJV) = sparse(ijv.i, ijv.j, ijv.v, ijv.n, ijv.n)

compress(ijv::IJV) = IJV(sparse(ijv))
