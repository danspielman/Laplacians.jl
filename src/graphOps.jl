#=

Operations on graphs.
Including transformations, combinations, and code for plotting.

Started by Dan Spielman.
Contributors: your name here.

Contains:

shortIntGraph - for converting the index type of a graph to an Int32.
lap - to produce the laplacian of a graph
unweight - change all the weights to 1
mapweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind},f) - to apply the function f to the weight of every edge.
uniformWeight : an example of mapweight.  It ignores the weight, and maps every weight to a random in [0,1]
productGraph(a0::SparseMatrixCSC, a1::SparseMatrixCSC)
edgeVertexMat(mat::SparseMatrixCSC) - signed edge vertex matrix
subsampleEdges(a::SparseMatrixCSC{Float64,Int64}, p::Float64)-
  produce a new graph that keeps each edge with probability p
twoLift(a, flip::Array)
function joinGraphs{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}, b::SparseMatrixCSC{Tval,Tind}, k::Integer)
plotGraph(gr,x,y,color=[0,0,1];dots=true,setaxis=true,number=false)
spectralDrawing(a)
toUnitVector(a::Array{Int64,1})
diagmat{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti})

=#


shortIntGraph(a::SparseMatrixCSC) = SparseMatrixCSC{Float64,Int32}(convert(Int32,a.m), convert(Int32,a.n), convert(Array{Int32,1},a.colptr), convert(Array{Int32,1},a.rowval), a.nzval)


lap(a) = spdiagm(a*ones(size(a)[1])) - a


function unweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind})
  (ai, aj, av) = findnz(a)
  n = size(a)[1]
  return sparse(ai,aj,ones(length(ai)),n,n)

end # unweight

# apply f to each nonzero entry of a
function mapweight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind},f)
    b = triu(a,1)
    b.nzval = map(f, b.nzval)
    b = b + b'

    return b

end # mapweight

# put a uniform [0,1] weight on every edge
uniformWeight{Tval,Tind}(a::SparseMatrixCSC{Tval,Tind}) = mapweight(a,x->rand(1)[1])

function productGraph(a0::SparseMatrixCSC, a1::SparseMatrixCSC)
  n0 = size(a0)[1]
  n1 = size(a1)[1]
  a = kron(speye(n0),a1) + kron(a0,speye(n1));

end # productGraph


# this produces the edge-vertex adjacency matrix
function edgeVertexMat(mat::SparseMatrixCSC)
    (ai,aj) = findnz(triu(mat,1))
    m = length(ai)
    n = size(mat)[1]
    return sparse(collect(1:m),ai,1.0,m,n) - sparse(collect(1:m),aj,1.0,m,n)
end



function subsampleEdges(a::SparseMatrixCSC{Float64,Int64}, p::Float64)
  (ai, aj, av) = findnz(a)
  n = size(a)[1]
  mask = map(x -> convert(Float64,x < p), rand(length(ai)))
  return sparse(ai,aj,mask.*av,n,n)

end # subsampleEdges


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
 create a disjoint union of graphs a and b,
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

function spectralDrawing(a)

    x, y = spectralCoords(a)
    plotGraph(a,x,y)

end # spectralDrawing


function spectralCoords(a)

    E = eigs(lap(a), nev = 3, which=:SR)
    V = E[2]
    return V[:,2], V[:,3]

end # spectralCoords

" creates a unit vector of length n from a given set of integers, with weights based on the number of occurences "
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

" returns the diagonal matrix(as a sparse matrix) of a graph "
function diagmat{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti})

  dw = zeros(Tv, G.n)

  u,v,w = findnz(G)

  for i in 1:length(u)
    dw[v[i]] = dw[v[i]] + w[i]
  end

  return sparse(collect(1:G.n), collect(1:G.n), dw)

end # diagmat

