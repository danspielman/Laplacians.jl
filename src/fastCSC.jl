#=

These are fast algorithms for maniuplating Julia sparse matrices.
It is part of the Laplacians package.

Started by Dan Spielman.

=#


""" Permute the rows and columns of a symmetric matrix, a little faster than the default
routine in Julia.  It basically does a radix sort on rows than columns. It does require `a`
to be symmetrix.  But, the same ideas could probably speed up many sparse matrix routines.
`sympermute(a, perm)` is the same as `a[perm,perm]`.
"""
function sympermute{Tv,Ti}(a::SparseMatrixCSC{Tv, Ti}, perm::Array{Ti,1})

    numnz = nnz(a)

    Ti1 = Ti
    Tv1 = Tv

    
    permi = Array(Ti1, a.n)
    permi[perm] = collect(1:a.n)
    
    J = Array(Ti1, numnz)

    deg = Array(Ti1, a.n)

    ptr = 1
    for col = 1:a.n
        deg[permi[col]] = a.colptr[col+1] - a.colptr[col]
        for k = a.colptr[col] : (a.colptr[col+1]-1)
            J[ptr] = permi[col]
            ptr += 1
        end
    end

    I = permi[a.rowval];

    
    I1 = Array(Ti1, numnz)
    J1 = Array(Ti1, numnz)
    V1 = Array(Tv1, numnz)

    cumdeg = cumsum(deg)
    col = [1;cumdeg+1]

    cumdeg1 = copy(cumdeg)
    
    for i in numnz:-1:1
        ptr = cumdeg1[I[i]]
        cumdeg1[I[i]] -= 1
        I1[ptr] = I[i]
        J1[ptr] = J[i]
        V1[ptr] = a.nzval[i]
    end

    I2 = Array(Ti1, numnz)
    V2 = Array(Tv1, numnz)

    for i in numnz:-1:1
        ptr = cumdeg[J1[i]]
        cumdeg[J1[i]] -= 1
        I2[ptr] = I1[i]
        V2[ptr] = V1[i]
    end

    aperm = SparseMatrixCSC(a.n, a.n, col, I2, V2)

    return aperm
end

#=
  This is a version of sympermute that is not quite as fast.
  But, it could be the basis of a faster algorithm.
  =#
function sympermuteB{Tv,Ti}(a::SparseMatrixCSC{Tv, Ti}, perm::Array{Ti,1})

    numnz = nnz(a)

    Ti1 = Ti
    Tv1 = Tv

    
    permi = Array(Ti1, a.n)
    permi[perm] = collect(1:a.n)
    
    J = Array(Ti1, numnz)
    I = Array(Ti1, numnz)
    deg = Array(Ti1, a.n)

    ptr = 1
    for col = 1:a.n
        deg[permi[col]] = a.colptr[col+1] - a.colptr[col]
        for k = a.colptr[col] : (a.colptr[col+1]-1)
            J[ptr] = permi[col]
            I[ptr] = permi[a.rowval[k]]
            ptr += 1
        end
    end
    
    I1 = Array(Ti1, numnz)
    J1 = Array(Ti1, numnz)
    V1 = Array(Tv1, numnz)

    cumdeg = cumsum(deg)
    for i in numnz:-1:1
        ptr = cumdeg[I[i]]
        I1[ptr] = I[i]
        J1[ptr] = J[i]
        V1[ptr] = a.nzval[i]
        cumdeg[I[i]] -= 1
    end

    I2 = Array(Ti1, numnz)
    V2 = Array(Tv1, numnz)
    cumdeg = cumsum(deg)
    for i in numnz:-1:1
        ptr = cumdeg[J1[i]]
        I2[ptr] = I1[i]
        V2[ptr] = V1[i]
        cumdeg[J1[i]] -= 1
    end

    cumdeg = cumsum(deg)
    col = [1;cumdeg+1]
    aperm = SparseMatrixCSC(a.n, a.n, col, I2, V2)

    return aperm
end


