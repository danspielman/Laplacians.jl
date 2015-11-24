

import Base.randperm


"The path graph on n vertices"
function pathGraph(n::Int64)
  x = append!(collect(1:(n-1)), collect(2:n))
  y = append!(collect(2:n), collect(1:(n-1)))
  w = ones(2 * (n - 1))
  return sparse(x, y, w)
end # pathGraph


"The complete graph" 
function completeGraph(n::Int64)
  return sparse(ones(n,n) - eye(n))
end # completeGraph


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

"""A random generalized ring graph of degree k.
Gens always contains 1, and the other k-1 edge types
are chosen from an exponential distribution"""
function randGenRing(n::Int64, k::Integer)
  gens = [1; 1 + ceil(Integer,exp(rand(k-1)*log(n-1)))]
  return generalizedRing(n, gens)
end


"""
~~~
generalizedNecklace{Tv, Ti}(A::SparseMatrixCSC{Tv, Ti}, H::SparseMatrixCSC, k::Int64)
~~~

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
    a = a + sparse(2:n,ceil(Integer,collect(1:n-1).*rand(n-1)),1,n,n)
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

        s = sort(ceil(Integer,n*rand(k)))
        good = (minimum(s[2:end]-s[1:(end-1)]) > 0)
        while good == false
            s = sort(ceil(Integer,n*rand(k)))
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


"""Randomly permutes the vertex indices"""
function randperm(mat::AbstractMatrix)
    perm = randperm(mat.n)
    return mat[perm,perm]
end

randperm(f::Expr) = randperm(eval(f))


"""Generate a random graph on n vertices with m edges.
The actual number of edges will probably be smaller, as we sample
with replacement"""
function ErdosRenyi(n::Integer, m::Integer)
    ai = rand(1:n, m)
    aj = rand(1:n, m)
    ind = (ai .!= aj)
    mat = sparse(ai[ind],aj[ind],1,n,n)
    uniformWeight!(mat)
    return mat
end

"""Generate an ER graph with average degree k,
and then return the largest component.
Will probably have fewer than n vertices.
If you want to add a tree to bring it back to n,
try ErdosRenyiClusterFix."""
function ErdosRenyiCluster(n::Integer, k::Integer)
    m = ceil(Integer,n*k/2)
    ai = rand(1:n, m)
    aj = rand(1:n, m)
    ind = (ai .!= aj)
    mat = sparse(ai[ind],aj[ind],1,n,n)
    mat = mat + mat'
   
    return biggestComp(mat)
end

"""Like an Erdos-Renyi cluster, but add back a tree so
it has n vertices"""
function ErdosRenyiClusterFix(n::Integer, k::Integer)
    m1 = ErdosRenyiCluster(n, k)
    n2 = n - size(m1)[1]
    if (n2 > 0)
        m2 = completeBinaryTree(n2)
        return joinGraphs(m1,m2,1)
    else
        return m1
    end
end

    
    
"""Generate a random graph with n vertices from one of our natural distributions"""
function pureRandomGraph(n::Integer)
   
    gr = []
    wt = []
    
    push!(gr,:(pathGraph($n)))
    push!(wt,1)

    push!(gr,:(ringGraph($n)))
    push!(wt,3)

    push!(gr,:(completeBinaryTree($n)))
    push!(wt,3)

    push!(gr,:(grownGraph($n,2)))
    push!(wt,6)

    push!(gr,:(grid2(ceil(Integer,sqrt($n)))[1:$n,1:$n]))
    push!(wt,6)

    push!(gr,:(randRegular($n,3)))
    push!(wt,6)

    push!(gr,:(ErdosRenyiClusterFix($n,2)))
    push!(wt,6)

    push!(gr,:(randGenRing($n,4)))
    push!(wt,6)

    i = sampleByWeight(wt)

    # make sure get a connected graph
    its = 0
    mat = eval(gr[i])

    while (~isConnected(mat)) && (its < 100)
        i = sampleByWeight(wt)
        mat = eval(gr[i])
        its += 1
    end
    if its == 100
        error("Getting a disconnected graph from $(gr[i])")
    end

    return floatGraph(mat)
      
end

"""sample an index with probability proportional to its weight given here"""
function sampleByWeight(wt)
    r = rand(1)*sum(wt)
    find(cumsum(wt) .> r)[1]
end

"""Builds a chimeric graph on n vertices.
The components come from pureRandomGraph,
connected by joinGraphs, productGraph and generalizedNecklace"""
function chimera(n::Integer)
    if (n < 30) || (rand() < .2)

        gr = pureRandomGraph(n)

        return gr
    end

    if (n < 200) 
        # just join disjoint copies of graphs

        n1 = 10 + floor(Integer,(n-20)*rand())
        n2 = n - n1
        k = ceil(Integer,exp(rand()*log(min(n1,n2)/2)))

        gr = joinGraphs(chimera(n1),chimera(n2),k)

        return gr
    end

    # split with probability .7

    if (rand() < .7)
        n1 = ceil(Integer,10*exp(rand()*log(n/20)))

        n2 = n - n1
        k = floor(Integer,1+exp(rand()*log(min(n1,n2)/2)))

        gr = joinGraphs(chimera(n1),chimera(n2),k)

        return gr

    else
        n1 = floor(Integer,10*exp(rand()*log(n/100)))

        n2 = floor(Integer, n / n1)

        if (rand() < .5)

            gr = productGraph(chimera(n1),chimera(n2))

        else

            k = floor(Integer,1+exp(rand()*log(min(n1,n2)/10)))
            gr = generalizedNecklace(chimera(n1),chimera(n2),k)


        end

        n3 = n - size(gr)[1]
        if (n3 > 0)
            gr = joinGraphs(gr,chimera(n3),2)
        end

        return gr
        
    end
end

"""Builds the kth chimeric graph on n vertices.
It does this by resetting the random number generator seed.
It should captute the state of the generator before that and then
return it, but it does not yet."""
function chimera(n::Integer, k::Integer)
    srand(k)
    g = chimera(n)
    return g
end

"""Applies one of a number of random weighting schemes to the edges of the graph"""
function randWeight(a)

    if (rand() < .1)
        return a
    end
    
    n = a.n
    (ai,aj) = findnz(a)
    m = length(ai)
    
    # potentials or edge-based

    if (rand() < .3)
        w = rand(m)

    else
        v = randn(a.n)

        # mult by matrix ?
        if (rand() < .5)

            invdeg = spdiagm(1./(a*ones(size(a)[1])))
            if (rand() < .5)
                for i in 1:10
                    v = a * (invdeg * v)
                    v = v - mean(v)
                end
            else
                for i in 1:10
                    v = v - a * (invdeg * v)
                    v = v - mean(v)
                end
            end
        end

        w = abs(v[ai]-v[aj]) 

    end

    # reciprocate or not?

    w[w.==0] = 1
    w[isnan(w)] = 1

    if (rand() < .5)
        w = 1./w
    end

    w = w / mean(w)

    ar = sparse(ai,aj,w,n,n)
    ar = ar + ar';
    return ar
end

"""Builds the kth wted chimeric graph on n vertices.
It does this by resetting the random number generator seed.
It should captute the state of the generator before that and then
return it, but it does not yet."""
function wtedChimera(n::Integer, k::Integer)
    srand(k)
    g = wtedChimera(n)
    return g
end

"""Generate a chimera, and then apply a random weighting scheme"""
function wtedChimera(n::Integer)
    return randWeight(chimera(n))
end

