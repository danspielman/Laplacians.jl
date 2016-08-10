
"Convert the indices in a graph to 32-bit ints.  This takes less storage, but does not speed up much"
shortIntGraph(a::SparseMatrixCSC) = SparseMatrixCSC{Float64,Int32}(convert(Int32,a.m), convert(Int32,a.n), convert(Array{Int32,1},a.colptr), convert(Array{Int32,1},a.rowval), a.nzval)

"Convert the nonzero entries in a graph to Float64"
floatGraph(a::SparseMatrixCSC) = SparseMatrixCSC{Float64,Int64}(a.m, a.n, a.colptr, a.rowval, convert(Array{Float64,1},a.nzval))



"""Create a Laplacian matrix from an adjacency matrix.
We might want to do this differently, say by enforcing symmetry"""
function lap{Tv,Ti}(a::SparseMatrixCSC{Tv,Ti})  
  return spdiagm(a * ones(Tv,size(a)[1])) - a
end

"""Create an adjacency matrix and a diagonal vector from a Laplacian with added diagonal weights"""
function adj{Tv,Ti}(la::SparseMatrixCSC{Tv,Ti})
  a = -(2 * la - triu(la) - tril(la))
  d = diag(la) - diag(lap(a))
  return a,d
end


"""Create a new graph in that is the same as the original, but with all edge weights 1"""
function unweight{Tval,Tind}(ain::SparseMatrixCSC{Tval,Tind})
    a = copy(ain)
    m = length(a.nzval)
    for i in 1:m
        a.nzval[i] = 1
    end
    return a
end # unweight

"""Change the weight of every edge in a to 1"""
function unweight!{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
    m = length(a.nzval)
    for i in 1:m
        a.nzval[i] = 1
    end
end # unweight


"""Create a new graph that is the same as the original, but with f applied to each nonzero entry of a. For example, to make the weight of every edge uniform in [0,1], we could write

~~~julia
b = mapweight(a, x->rand(1)[1])
~~~
"""
function mapweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind},f)
    b = triu(a,1)
    b.nzval = map(f, b.nzval)
    b = b + b'

    return b

end # mapweight

"""Put a uniform [0,1] weight on every edge.  This is an example of how to use mapweight."""
uniformWeight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}) = mapweight(a,x->rand(1)[1])

"""Set the weight of every edge to 1"""
function uniformWeight!(mat::SparseMatrixCSC)
    mat.nzval = ones(length(mat.nzval))
end

"""The Cartesian product of two graphs.  When applied to two paths, it gives a grid."""
function productGraph(a0::SparseMatrixCSC, a1::SparseMatrixCSC)
  n0 = size(a0)[1]
  n1 = size(a1)[1]
  a = kron(speye(n0),a1) + kron(a0,speye(n1));

end # productGraph


"""The signed edge-vertex adjacency matrix"""
function edgeVertexMat(mat::SparseMatrixCSC)
    (ai,aj) = findnz(triu(mat,1))
    m = length(ai)
    n = size(mat)[1]
    return sparse(collect(1:m),ai,1.0,m,n) - sparse(collect(1:m),aj,1.0,m,n)
end


"""Create a new graph from the old, but keeping edge edge with probability `p`"""
function subsampleEdges(a::SparseMatrixCSC{Float64,Int64}, p::Float64)
  (ai, aj, av) = findnz(triu(a))
  n = size(a)[1]
  mask = map(x -> convert(Float64,x < p), rand(length(ai)))
  as = sparse(ai,aj,mask.*av,n,n)
  as = as + as'
  return as

end # subsampleEdges


"""Creats a 2-lift of a.  `flip` is a boolean indicating which edges cross"""
function twoLift(a, flip::AbstractArray{Bool,1})
    (ai,aj,av) = findnz(triu(a))
    m = length(ai)
    #flip = rand(false:true,m)
    n = size(a)[1]
    a0 = sparse(ai[flip],aj[flip],1,n,n)
    a1 = sparse(ai[!flip],aj[!flip],1,n,n)
    a00 = a0 + a0'
    a11 = a1 + a1'
    return [a00 a11; a11 a00]
end

twoLift(a) = twoLift(a,rand(false:true,div(nnz(a),2)))

twoLift(a, k::Integer) = twoLift(a,randperm(div(nnz(a),2)) .> k)




"""
 Create a disjoint union of graphs a and b,
 and then put k random edges between them
"""
function joinGraphs{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, b::SparseMatrixCSC{Tval,Tind}, k::Integer)
    na = size(a)[1]
    nb = size(b)[1]

    (ai,aj,av) = findnz(a)
    (bi,bj,bv) = findnz(b)
    bi = bi + na
    bj = bj + na

    ji = rand(1:na,k)
    jj = rand(1:nb,k)+na

    ab = sparse([ai;bi;ji;jj],[aj;bj;jj;ji],[av;bv;ones(Tval,2*k)],na+nb,na+nb)
end


"""
 Create a disjoint union of graphs a and b,
  with no edges between them.
"""
disjoin(a,b) = joinGraphs(a,b,0)



"""Plots graph gr with coordinates (x,y)"""
function plotGraph(gr,x,y,color=[0,0,1];dots=true,setaxis=true,number=false)
  (ai,aj,av) = findnz(triu(gr))
  arx = [x[ai]';x[aj]';NaN*ones(length(ai))']
  ary = [y[ai]';y[aj]';NaN*ones(length(ai))']
  p = plot(arx[:],ary[:],color=color)

  if dots
    plot(x,y,color=color,marker="o",linestyle="none")
  end #if


  if number
    for i in 1:length(x)
      annotate(i, xy = [x[i]; y[i]])
    end
  end

  axis("off")

  if setaxis
  ax = axes()

  minx = minimum(x)
  maxx = maximum(x)
  miny = minimum(y)
  maxy = maximum(y)
  delx = maxx - minx
  dely = maxy - miny
  ax[:set_ylim]([miny - dely/20, maxy + dely/20])
  ax[:set_xlim]([minx - delx/20, maxx + delx/20])
  end

  return p
end # plotGraph

"""Computes spectral coordinates, and then uses plotGraph to draw"""
function spectralDrawing(a)

    x, y = spectralCoords(a)
    plotGraph(a,x,y)

end # spectralDrawing

"""Computes the spectral coordinates of a graph"""
function spectralCoords(a)

    E = eigs(lap(a), nev = 3, which=:SR)
    V = E[2]
    return V[:,2], V[:,3]

end # spectralCoords

"""Creates a unit vector of length n from a given set of integers, with weights based on the number of occurences"""
function toUnitVector(a::Array{Int64,1}, n)

  v = zeros(n)

  for i in a
    v[i] = v[i] + 1
  end

  sum = 0
  for i in 1:length(v)
    sum = sum + v[i] * v[i]
  end

  for i in 1:length(v)
    v[i] = v[i] / sqrt(sum)
  end

  return v

end # toUnitVector

"Returns the diagonal matrix(as a sparse matrix) of a graph"
function diagmat{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti})

  dw = zeros(Tv, G.n)

  u,v,w = findnz(G)

  for i in 1:length(u)
    dw[v[i]] = dw[v[i]] + w[i]
  end

  return sparse(collect(1:G.n), collect(1:G.n), dw)

end # diagmat


"""
 Constructs a generalized necklace graph starting with two graphs A and H. The
resulting new graph will be constructed by expanding each vertex in H to an
instance of A. k random edges will be generated between components. Thus, the
resulting graph may have weighted edges.
"""
function generalizedNecklace{Tv, Ti}(A::SparseMatrixCSC{Tv, Ti}, H::SparseMatrixCSC, k::Int64)
  a = findnz(A)
  h = findnz(H)

  # these are square matrices
  n = A.n
  m = H.n

  newI = Ti[]
  newJ = Ti[]
  newW = Tv[]

  # duplicate the vertices in A so that each vertex in H corresponds to a copy of A
  for i in 1:m
    newI = append!(newI, a[1] + n * (i - 1))
    newJ = append!(newJ, a[2] + n * (i - 1))
    newW = append!(newW, a[3])
  end

  # for each edge in H, add k random edges between two corresponding components
  # multiedges will be concatenated to a single edge with higher cost
  for i in 1:length(h[1])
    u = h[1][i]
    v = h[2][i]

    if (u < v)
      #component x is from 1 + (x - 1) * n to n + (x - 1) * n
      for edgeToAdd in 1:k
        newU = rand(1:n) + n * (u - 1)
        newV = rand(1:n) + n * (v - 1)
        append!(newI, [newU, newV])
        append!(newJ, [newV, newU])
        append!(newW, [1, 1])
      end
    end
  end

  return sparse(newI, newJ, newW)
end # generalizedNecklace
