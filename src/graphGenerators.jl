
#=

Started by Dan Spielman.
Contributors: ???

Code for reading graphs from files, and for constructing special graphs.

Right now, it has:

 readIJ(filename::AbstractString)
 readIJV(filename::AbstractString)
 writeIJV(filename::AbstractString, mat)
 ringGraph(n::Int64)
 generalizedRing(n::Int64, gens)
 randMatching(n::Int64)
 randRegular(n::Int64, k::Int64)
 grownGraph(n::Int64, k::Int64)
 grownGraphD(n::Int64, k::Int64)
 prefAttach(n::Int64, k::Int64, p::Float64)
 hyperCube(d::Int64)
 completeBinaryTree(n::Int64)
 grid2(n::Int64)
 grid2(n::Int64, m::Int64; isotropy=1)
 grid2coords(n::Int64, m::Int64)
 generalizedNecklace(A::SparseMatrixCSC, H::SparseMatrixCSC, k::Int64)
 completeGraph(n::Int64)
 pathGraph(n::Int64)
=#


# to read a simple edge list, each line being an (i, j) pair
function readIJ(filename::AbstractString)
  edlist = readdlm(filename,',')
  n = maximum(edlist)
  m = size(edlist)
  edlist = convert(Array{Int64,2}, edlist)
  a = sparse(edlist[:,1],edlist[:,2],ones(m[1]),n,n)
  a = a + a'

  return a
end # readIJ

function readIJV(filename::AbstractString)
  data = readdlm(filename,',')
  n = maximum(data[:,1:2])
  m = size(data)
  edlist = convert(Array{Int64,2}, data[:,1:2])
  wts = convert(Array{Float64,1}, data[:,3])

  a = sparse(edlist[:,1],edlist[:,2],wts,n,n)
  a = a + a'

  return a
end # readIJV


# just writes the upper triangular portion of it
function writeIJV(filename::AbstractString, mat)

  (ai,aj,av) = findnz(triu(mat))
  fh = open(filename,"w")
  for i in 1:length(ai)
    write(fh, "$(ai[i]),$(aj[i]),$(av[i])\n")
  end
  close(fh)

end #write IJV



"""The simple ring on n vertices"""
function ringGraph(n::Int64)
    a = spdiagm(ones(n-1),1,n,n)
    a[1,n] = 1
    a = a + a'
end

"""A generalization of a ring graph.
The vertices are integers modulo n.
Two are connected if their difference is in gens.
For example, 

```
generalizedRing(17, [1 5])
```
"""
function generalizedRing(n::Int64, gens)
    k = length(gens)
    m = 2*n*k
    ai = zeros(Int64,m)
    aj = zeros(Int64,m)
    ind = 1
    for i in 0:(n-1)
        for j in 1:k
            ai[ind] = i
            aj[ind] = mod(i+gens[j],n)
            ind = ind + 1
            ai[ind] = i
            aj[ind] = mod(i-gens[j],n)
            ind = ind + 1
        end
    end
    return sparse(1+ai,1+aj,ones(m),n,n)
    #return ai, aj
end


"""A random matching on n vertices"""
function randMatching(n::Int64)

  p = randperm(n)
  n1 = convert(Int64,floor(n/2))
  n2 = 2*n1
  a = sparse(p[1:n1],p[(n1+1):n2],ones(n1),n,n)

  a = a + a'

  return a

end # randMatching

"""A sum of k random matchings on n vertices"""
function randRegular(n::Int64, k::Int64)
  a = randMatching(n)
  for i in 2:k
    a = a + randMatching(n)
  end

  return a
end # randRegular


"""Create a graph on n vertices.
For each vertex, give it k edges to randomly chosen prior
vertices.
This is a variety of a preferential attachment graph.    
"""
function grownGraph(n::Int64, k::Int64)
  a = spzeros(n,n)

  for i = 1:k
    a = a + sparse(2:n,iceil(collect(1:n-1).*rand(n-1)),1,n,n)
  end

  a = a + a'
end # grownGraph

# used in grownGraphD
function randSet(n::Integer,k::Integer)
    if n == k
        return collect(1:n)
    elseif n < k
        error("n must be at least k")
    else

        s = sort(iceil(n*rand(k)))
        good = (minimum(s[2:end]-s[1:(end-1)]) > 0)
        while good == false
            s = sort(iceil(n*rand(k)))
            good = (minimum(s[2:end]-s[1:(end-1)]) > 0)
        end

        return s

    end
end

"""Like a grownGraph, but it forces the edges to all be distinct.
It starts out with a k+1 clique on the first k vertices"""
function grownGraphD(n::Int64, k::Int64)
    a = spzeros(n,n)

    u = zeros(Int64, k*(n-k-1))
    v = zeros(Int64, k*(n-k-1))

    for i in (k+2):n
        nb = randSet(i-1,k)
        u[(i-k-2)*k + collect(1:k)] = i
        v[(i-k-2)*k + collect(1:k)] = nb
    end

    a = sparse(u,v,1,n,n)

    (ai,aj) = findnz(triu(ones(k+1,k+1),1))
    a = a + sparse(ai,aj,1,n,n)
    a = a + a'

end # grownGraphD

"""A preferential attachment graph in which each vertex has k edges to those
that come before.  These are chosen with probability p to be from a random vertex,
and with probability 1-p to come from the endpoint of a random edge.
It begins with a k-clique on the first k+1 vertices."""
function prefAttach(n::Int64, k::Int64, p::Float64)
    if n == (k+1)
        return sparse(ones(Float64,n,n) - eye(Float64,n))
    elseif n <= k
        error("n must be more than k")
    else

        u = zeros(Int64,n*k)
        v = zeros(Int64,n*k)


        # fill in the initial clique
        # this will accidentally double every edge in the clique
        # we clean it up at the end
        ind = 1
        for i in 1:(k+1)
            for j in 1:(k+1)
                if i != j
                    u[ind] = i
                    v[ind] = j
                    ind += 1
                end
            end
        end

        s = zeros(Int64,k)
        for i in (k+2):n
            distinct = false
            while distinct == false
                for j in 1:k
                    if rand(Float64) < p
                        s[j] = rand(1:(i-1))
                    else
                        s[j] = v[rand(1:(k*(i-1)))]
                    end
                end
                s = sort(s)
                distinct = true
                for ii in 1:(k-1)
                    if s[ii] == s[ii+1]
                        distinct = false
                    end
                end
                # distinct = (minimum(s[2:end]-s[1:(end-1)]) > 0)

            end

            for j in 1:k
                u[ind] = i
                v[ind] = s[j]
                ind += 1
            end

        end # for i

        w = ones(Float64,n*k)

        w[1:(k*(k+1))] = 1/2


        a = sparse(u,v,w,n,n)
        a = a + a'
        return a
    end
end

"""The d dimensional hypercube.  Has 2^d vertices"""
function hyperCube(d::Int64)
  a = sparse([0 1; 1 0])

  for i = 1:(d-1)
    k = 2^i
    D = speye(k)
    a = [a D; D a]
  end

  return a
end # hyperCube

"""The complete binary tree on n vertices"""
function completeBinaryTree(n::Int64)
  k = div(n-1,2)
  a = sparse(collect(1:k),2*collect(1:k),1,n,n) + sparse(collect(1:k),2*collect(1:k)+1,1,n,n)

  if 2*k+1 < n
    a[n-1,n] = 1
  end

  a = a + a'

  return a
end # completeBinaryTree

"""An n-by-m grid graph.  iostropy is the weighting on edges in one direction."""
function grid2(n::Int64, m::Int64; isotropy=1)
  a = kron(speye(n),spdiagm(ones(m-1),1,m,m))
  a = a + isotropy*kron(spdiagm(ones(n-1),1,n,n), speye(m))
  a = a + a'
  return a
end # grid2

grid2(n::Int64) = grid2(n,n)

"""Coordinates for plotting the vertices of the n-by-m grid graph"""
function grid2coords(n::Int64, m::Int64)
  x = kron(collect(1:n),ones(m))
  y = kron(ones(n),collect(1:m))
  return x, y
end # grid2coords

grid2coords(n) = grid2coords(n, n)


" The complete graph " 
function completeGraph(n::Int64)
  return sparse(ones(n,n) - eye(n))
end # completeGraph

" The path graph on n vertices "
function pathGraph(n::Int64)
  x = append!(collect(1:(n-1)), collect(2:n))
  y = append!(collect(2:n), collect(1:(n-1)))
  w = ones(2 * (n - 1))
  return sparse(x, y, w)
end # pathGraph
