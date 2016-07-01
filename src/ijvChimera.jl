type IJV{Tv,Ti}
    n::Ti
    nnz::Ti
    i::Array{Ti,1}
    j::Array{Ti,1}
    v::Array{Tv,1}
end


function IJV{Tv,Ti}(A::SparseMatrixCSC{Tv, Ti}) 
    (ai,aj,av) = findnz(A)
    IJV{Tv,Ti}(A.n, nnz(A), ai, aj, av)
end


import Base.SparseMatrix.sparse
sparse(IJV) = sparse(IJV.i, IJV.j, IJV.v, IJV.n, IJV.n)
sparse(ijv)
sum(abs(a-sparse(ijv)))

function generalizedNecklace{Tv, Ti}(A::IJV{Tv, Ti}, H::IJV{Tv, Ti}, k::Int64)

  # these are square matrices
  n = A.n
  m = H.n

  newI = Ti[]
  newJ = Ti[]
  newW = Tv[]

  # duplicate the vertices in A so that each vertex in H corresponds to a copy of A
  for i in 1:m
    newI = append!(newI, A.i + n * (i - 1))
    newJ = append!(newJ, A.j + n * (i - 1))
    newW = append!(newW, A.v)
  end

  # for each edge in H, add k random edges between two corresponding components
  # multiedges will be concatenated to a single edge with higher cost
  for i in 1:H.nnz
    u = H.i[i]
    v = H.j[i]

    if (u < v)
      #component x is from 1 + (x - 1) * n to n + (x - 1) * n
      for edgeToAdd in 1:k
                r1 = rand(1:n)
                r2 = rand(1:n)
        newU = r1 + n * (u - 1)
        newV = r2 + n * (v - 1)
        append!(newI, [newU, newV])
        append!(newJ, [newV, newU])
        append!(newW, [1, 1])
      end
    end
  end

    return IJV{Tv, Ti}(n*m, length(newI), newI, newJ, newW)
end # generalizedNecklace


function generalizedNecklace2{Tv, Ti}(A::IJV{Tv, Ti}, H::IJV{Tv, Ti}, k::Int64)

    # these are square matrices
    n = A.n
    m = H.n

    newI = Array(Ti, m*A.nnz)
    newJ = Array(Ti, m*A.nnz)
    newW = Array(Tv, m*A.nnz)
    
    ptr = 1
    for i in 1:m
        for j in 1:A.nnz
            newI[ptr] = A.i[j] + n*(i-1)
            newJ[ptr] = A.j[j] + n*(i-1)
            newW[ptr] = A.v[j] 
            ptr += 1
        end
    end

    

  # for each edge in H, add k random edges between two corresponding components
  # multiedges will be concatenated to a single edge with higher cost
    
    # this is not for efficiency, but for backwards compatability
    rnd = rand(1:n, k*H.nnz)
    nextI = rnd[1:2:end]
    nextJ = rnd[2:2:end]
    
    nextW = ones(Tv, k*H.nnz)
    
    ptr = 1
    for i in 1:H.nnz
        u = H.i[i]
        v = H.j[i]

        if (u < v)
            #component x is from 1 + (x - 1) * n to n + (x - 1) * n
            for edgeToAdd in 1:k
                nextI[ptr] += n * (u - 1)
                nextJ[ptr] += n * (v - 1)
                ptr += 1
            end
        end
    end
  
    newI = [newI; nextI; nextJ]
    newJ = [newJ; nextJ; nextI]
    newW = [newW; nextW]
    
    return IJV{Tv, Ti}(n*m, length(newI), newI, newJ, newW)
end # generalizedNecklace


