
"""
    graph = shortIntGraph(a::SparseMatrixCSC)

Convert the indices in a graph to 32-bit ints.
This takes less storage, but does not speed up much.
"""
shortIntGraph(a::SparseMatrixCSC) = SparseMatrixCSC{Float64,Int32}(convert(Int32,a.m), convert(Int32,a.n), convert(Array{Int32,1},a.colptr), convert(Array{Int32,1},a.rowval), a.nzval)


"""
  graph = floatGraph(a::SparseMatrixCSC)

Convert the nonzero entries in a graph to Float64.
"""
floatGraph(a::SparseMatrixCSC) = SparseMatrixCSC{Float64,Int64}(a.m, a.n, a.colptr, a.rowval, convert(Array{Float64,1},a.nzval))



"""
    l = lap(a)

Create a Laplacian matrix from an adjacency matrix. We might want to do this differently, say by enforcing symmetry
"""
lap(a) = sparse(Diagonal(a*ones(size(a)[1]))) - a

"""
    a,d = adj(sddm)

Create an adjacency matrix and a diagonal vector from an SDD M-matrix.
That is, from a Laplacian with added diagonal weights
"""
function adj(sddm::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    a = sparse(Diagonal(diag(sddm))) - sddm
    d = sddm*ones(size(sddm,1))
    return a,d
end

function adj(sddm::Array{Tv,2}) where Tv
    return adj(sparse(sddm))
end


"""
Add a new vertex to a with weights to the other vertices corresponding to diagonal surplus weight.

This is an efficient way of writing [a d; d' 0]
"""
function extendMatrix(a::SparseMatrixCSC{Tv,Ti}, d::Array{Tv,1}) where {Tv,Ti}

    @assert size(a,1) == length(d)

    if sum(abs.(d)) == 0
        return a
    end
    
    dpos = d.*(d.>0)

    n = length(d)
    
    ai,aj,av=findnz(a)
    
    ai2 = [ai;1:n;(n+1)*ones(Int,n)]
    aj2 = [aj;(n+1)*ones(Int,n);1:n]
    av2 = [av;dpos;dpos]
    a2 = sparse(ai2,aj2,av2)
    
    return a2
    
end


"""
    wt1 = unweight(a)

Create a new graph in that is the same as the original, but with all edge weights 1
"""
function unweight(ain::SparseMatrixCSC{Tval,Tind}) where {Tval,Tind}
    a = copy(ain)
    m = length(a.nzval)
    for i in 1:m
        a.nzval[i] = one(Tval)
    end
    return a
end # unweight

"""
    unweight!(a)

Change the weight of every edge in a to 1
"""
function unweight!(a::SparseMatrixCSC{Tval,Tind}) where {Tval,Tind}
    m = length(a.nzval)
    for i in 1:m
        a.nzval[i] = one(Tval)
    end
end # unweight


function unweight!(ijv::IJV)
    ijv.v .= 1.0
end


"""
    b = mapweight(a, x->rand())

Create a new graph that is the same as the original, but with f applied to each nonzero entry of a. For example, to make the weight of every edge uniform in [0,1], we could write.
"""
function mapweight(a::SparseMatrixCSC{Tval,Tind},f) where {Tval,Tind}
    b = triu(a,1)
    for i in 1:length(b.nzval)
        b.nzval[i] = f(b.nzval[i])
    end
    b = b + b'

    return b

end # mapweight

"""
    wted = uniformWeight(unwted)

Put a uniform [0,1] weight on every edge.  This is an example of how to use mapweight."""
uniformWeight_ver(::Type{V06}, a::SparseMatrixCSC)  = mapweight(a,x->rand_ver(V06))
uniformWeight_ver(::Type{Vcur}, a::SparseMatrixCSC)  = mapweight(a,x->rand_ver(Vcur))
uniformWeight(a::SparseMatrixCSC)  = uniformWeight_ver(Vcur, a::SparseMatrixCSC)


"""
    uniformWeight!(a)

Set the weight of every edge to random uniform [0,1]
"""
uniformWeight!(mat::SparseMatrixCSC) = uniformWeight_ver!(Vcur, mat)

function uniformWeight_ver!(ver, mat::SparseMatrixCSC)
    for i in 1:length(mat.nzval)
        mat.nzval[i] = rand_ver(ver)
    end
end


"""
    aprod = productGraph(a0, a1)

The Cartesian product of two graphs.  When applied to two paths, it gives a grid.
"""
function product_graph(a0::SparseMatrixCSC, a1::SparseMatrixCSC)
  n0 = size(a0)[1]
  n1 = size(a1)[1]
  a = kron(sparse(I,n0,n0),a1) + kron(a0,sparse(I, n1, n1));

end # productGraph

function product_graph(b::IJV{Tva,Tia}, a::IJV{Tvb,Tib}) where {Tva, Tvb, Tia, Tib}

    Ti = promote_type(Tia, Tib)

    n = a.n * b.n

    @assert length(a.i) == a.nnz

    a_edge_from = kron(ones(Ti, a.nnz), a.n*collect(0:(b.n-1)))
    ai = a_edge_from + kron(a.i, ones(Ti, b.n))
    aj = a_edge_from + kron(a.j, ones(Ti, b.n))
    av = kron(a.v, ones(b.n))

    b_edge_from = kron(collect(1:a.n), ones(Ti, b.nnz))
    bi = b_edge_from + kron(ones(Ti, a.n), (b.i .- 1) .* a.n)
    bj = b_edge_from + kron(ones(Ti, a.n), (b.j .- 1) .* a.n)
    bv = kron(ones(a.n), b.v)

    return IJV(n, length(av)+length(bv),
        [ai; bi], [aj; bj], [av; bv])
end


"""
    ap = power(a::SparseMatrixCSC, k::Int)

Returns the kth power of a.
"""
function power(a::SparseMatrixCSC, k::Int)
  ap = a^k
  ap = ap - sparse(Diagonal(diag(ap)))
end


"""
    a_new = thicken_once(a)

Creates one edge for every vertex in a of degree > 1
by connecting two of its random neighbors.
To use this to thicken a, return unweight(a + a_new).

```
a = grid2(5)
a2 = unweight(a + thicken_once(a))
(x,y) = grid2coords(5,5);
plotGraph(a2,x,y)
```
"""
function thicken_once(a::SparseMatrixCSC; ver = Vcur)
    n = a.n
    e_new = zeros(Int,n,2)
    ptr = 0
    for i in 1:n
        d = deg(a, i)
        if d > 1
            nbr1 = 0
            nbr2 = 0
            while nbr1 == nbr2
                nbr1 = rand_ver(ver, 1:d)
                nbr2 = rand_ver(ver, 1:d)
            end
            ptr += 1
            e_new[ptr,:] = [nbri(a,i,nbr1) nbri(a,i,nbr2)]
        end
    end

    a_new = sparse(e_new[1:ptr,1], e_new[1:ptr,2], 1, n, n)
    return unweight(a_new + a_new')
end


"""
    a_new = thicken(A,k)

Create a new graph with at least k times as many edges as A
By connecting nodes with common neighbors at random.
When this stops working (not enough new edges),
repeat on the most recently produced graph.
If k is too big, it is decreased so the average degree will not
be pushed much above n/2.


When called without k, it just runs thicken_once.

For example:
```
a = grid2(5)
a2 = thicken(a,3)
(x,y) = grid2coords(5,5);
plotGraph(a2,x,y)
```
"""
function thicken(a::SparseMatrixCSC,k; ver = Vcur)
    n = a.n

    k = min(k, round(Int, n^2/nnz(a)/2 ))

    ne0 = nnz(a)/2

    ne = ne0

    a_new = a

    while nnz(a_new)/2 < ne0*k

        before = nnz(a_new)/2

        a_new = unweight(a_new + thicken_once(a, ver=ver))
        if nnz(a_new)/2 - before < n/4
            a = a_new
        end
    end
    return a_new

end

thicken(a; ver=Vcur) = unweight(a + thicken_once(a, ver=ver))

"""
    graph = generalizedNecklace(A, H, k::Int64)

Constructs a generalized necklace graph starting with two graphs A and H. The
resulting new graph will be constructed by expanding each vertex in H to an
instance of A. k random edges will be generated between components. Thus, the
resulting graph may have weighted edges.
"""
function generalizedNecklace(A::SparseMatrixCSC{Tv, Ti}, H::SparseMatrixCSC, k::Int64; ver=Vcur) where {Tv, Ti}
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
    newI = append!(newI, a[1] .+ n * (i - 1))
    newJ = append!(newJ, a[2] .+ n * (i - 1))
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
        newU = rand_ver(ver,1:n) + n * (u - 1)
        newV = rand_ver(ver,1:n) + n * (v - 1)
        append!(newI, [newU, newV])
        append!(newJ, [newV, newU])
        append!(newW, [1, 1])
      end
    end
  end

  return sparse(newI, newJ, newW)
end # generalizedNecklace

function generalized_necklace(A::IJV, H::IJV, k::Integer; ver=Vcur) 
  
    if ver == V06
        H = IJV(sparse(H))
    end

    # these are square matrices
    n = A.n
    m = H.n
  
    # duplicate the vertices in A so that each vertex in H corresponds to a copy of A
    newW = kron(ones(m), A.v)
    newI = kron(ones(Int, m), A.i) + kron(n * collect(0:(m-1)), ones(Int, A.nnz) )
    newJ = kron(ones(Int, m), A.j) + kron(n * collect(0:(m-1)), ones(Int, A.nnz) )
  
    # for each edge in H, add k random edges between two corresponding components
    # multiedges will be concatenated to a single edge with higher weight

    Hi = Array{Int}(undef, div(k*H.nnz,2))
    Hj = Array{Int}(undef, div(k*H.nnz,2))
    Hv = ones(div(k*H.nnz,2))

    ind = 0
    for i in 1:H.nnz
      u = H.i[i]
      v = H.j[i] 
  
      if (u < v)
        #component x is from 1 + (x - 1) * n to n + (x - 1) * n
        for j in 1:k
          newU = rand_ver(ver, 1:n) + n * (u - 1)
          newV = rand_ver(ver, 1:n) + n * (v - 1)

          ind += 1
          Hi[ind] = newU
          Hj[ind] = newV

        end
      end
    end
  
    return IJV(n*m, length(newI) + k*H.nnz, [newI; Hi; Hj] , [newJ; Hj; Hi], [newW; Hv; Hv])

  end # generalizedNecklace
  


"""
    U = edgeVertexMat(a)

The signed edge-vertex adjacency matrix
"""
function edgeVertexMat(mat::SparseMatrixCSC)
    (ai,aj) = findnz(triu(mat,1))
    m = length(ai)
    n = size(mat)[1]
    return sparse(collect(1:m),ai,1.0,m,n) - sparse(collect(1:m),aj,1.0,m,n)
end

"""
    U = wtedEdgeVertexMat(a)

The signed and weighted edge-vertex adjacency matrix, so U'*U = L
"""
function wtedEdgeVertexMat(mat::SparseMatrixCSC)
    (ai,aj,av) = findnz(triu(mat,1))
    m = length(ai)
    n = size(mat)[1]
    v = av.^(1/2)
    return sparse(collect(1:m),ai,v,m,n) - sparse(collect(1:m),aj,v,m,n)
end

"""
    graph = subsampleEdges(a::SparseMatrixCSC, p::Float64)

Create a new graph from the old, but keeping edge edge with probability `p`
"""
function subsampleEdges(a::SparseMatrixCSC, p::Float64; ver=Vcur)
  (ai, aj, av) = findnz(triu(a))
  n = size(a)[1]
  mask = rand_ver(ver, length(ai)) .< p
  as = sparse(ai[mask],aj[mask],av[mask],n,n)
  as = as + as'
  return as

end # subsampleEdges


"""
    graph = twoLift(a, flip::AbstractArray{Bool,1})
    graph = twoLift(a)
    graph = twoLift(a, k::Integer)

Creats a 2-lift of a.  `flip` is a boolean indicating which edges cross.
In the third version, k is the number of edges that cross.
"""
function twoLift(a, flip::AbstractArray{Bool,1})
    (ai,aj,av) = findnz(triu(a))
    m = length(ai)
    n = size(a)[1]
    a0 = sparse(ai[flip],aj[flip],1,n,n)
    a1 = sparse(ai[.!(flip)],aj[.!(flip)],1,n,n)
    a00 = a0 + a0'
    a11 = a1 + a1'
    return [a00 a11; a11 a00]
end

twoLift(a; ver=Vcur) = twoLift(a,rand_ver(ver, false:true,div(nnz(a),2)))

twoLift(a, k::Integer; ver=Vcur) = twoLift(a,randperm_ver(ver, div(nnz(a),2)) .> k)




"""
    graph = joinGraphs(a, b, k::Integer)

 Create a disjoint union of graphs a and b,
 and then put k random edges between them
"""
function join_graphs(a::SparseMatrixCSC{Tval,Tind}, b::SparseMatrixCSC{Tval,Tind}, k::Integer; ver=Vcur) where {Tval,Tind}
    na = size(a)[1]
    nb = size(b)[1]

    (ai,aj,av) = findnz(a)
    (bi,bj,bv) = findnz(b)
    bi = bi .+ na
    bj = bj .+ na

    ji = rand_ver(ver, 1:na,k)
    jj = rand_ver(ver, 1:nb,k) .+ na

    ab = sparse([ai;bi;ji;jj],[aj;bj;jj;ji],[av;bv;ones(Tval,2*k)],na+nb,na+nb)
end


function join_graphs(a::IJV, b::IJV, k::Integer; ver=Vcur) 

    ji = rand_ver(ver, 1:a.n,k)
    jj = rand_ver(ver, 1:b.n,k) .+ a.n

    gi = [a.i; b.i .+ a.n ; ji ; jj] 
    gj = [a.j; b.j .+ a.n ; jj; ji]
    gv = [a.v ; b.v ; ones(2*k)]

    return IJV(a.n + b.n, a.nnz + b.nnz + 2*k, gi, gj, gv)
end

"""
    graph = join_graphs!(a::IJV, b::IJV, k::Integer)

 Create a disjoint union of graphs a and b,
 and then put k random edges between them, merging b into a.
"""
function join_graphs!(a::IJV, b::IJV, k::Integer; ver=Vcur) 

    ji = rand_ver(ver, 1:a.n,k)
    jj = rand_ver(ver, 1:b.n,k) .+ a.n

    append!(a.i, [b.i .+ a.n ; ji ; jj] )
    append!(a.j, [b.j .+ a.n ; jj; ji])
    append!(a.v, [b.v ; ones(2*k)])
    a.n += b.n
    a.nnz += b.nnz + 2*k

end

"""
    graph = disjoin(a,b)

 Create a disjoint union of graphs a and b,
  with no edges between them.
"""
disjoin(a::SparseMatrixCSC,b::SparseMatrixCSC) = join_graphs(a,b,0)

disjoin(a::IJV, b::IJV) = IJV(a.n+b.n, a.nnz + b.nnz, 
        [a.i ; b.i .+ a.n],
        [a.j ; b.j .+ a.n],
        [a.v;b.v])

function disjoin!(a::IJV, b::IJV) 
    append!(a.i, b.i .+ a.n )
    append!(a.j, b.j .+ a.n )
    append!(a.v, b.v)
    a.n += b.n
    a.nnz += b.nnz
end

"""
    b = firstn(a::IJV, n::Integer)

Only keep the first n vertices of a.
"""
function firstn(a::IJV, n::Integer) 
    mask = (a.i .<= n) .& (a.j .<= n)
    return IJV(n, sum(mask), a.i[mask], a.j[mask], a.v[mask])
end





"""
    plotGraph(gr,x,y,color=[0,0,1];dots=true,setaxis=true,number=false)

Plots graph gr with coordinates (x,y)
"""
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

"""
    spectralDrawing(a)

Computes spectral coordinates, and then uses plotGraph to draw
"""
function spectralDrawing(a)

    x, y = spectralCoords(a)
    plotGraph(a,x,y)

end # spectralDrawing

"""
    spectralCoords(a)

Computes the spectral coordinates of a graph
"""
function spectralCoords(a)

    E = eigs(lap(a), nev = 3, which=:SR)
    V = E[2]
    return V[:,2], V[:,3]

end # spectralCoords


"""
    d = diagmat(a)

Returns the diagonal weighted degree matrix(as a sparse matrix) of a graph
"""
function diagmat(a::SparseMatrixCSC{Tv, Ti}) where {Tv, Ti}

  return sparse(Diagonal(vec(sum(a,dims=1))))

end # diagmat
