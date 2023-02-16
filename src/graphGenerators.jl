

import Random.randperm
import RandomV06.randperm_ver


"""
    ijv = empty_graph_ijv(n)
"""
empty_graph_ijv(n::Integer) = IJV{Float64,Int64}(Int64(n),0,[],[],[])

"""
    ijv = empty_graph(n)
"""
empty_graph(n::Integer) = sparse(empty_graph_ijv(n))


"""
    graph = path_graph(n)
"""
path_graph(n) = sparse(path_graph_ijv(n))

"""
    ijv = path_graph_ijv(n::Int64)
"""
function path_graph_ijv(n::Integer)
    IJV(n, 2*(n-1),
        [collect(1:(n-1)) ; collect(2:n)], 
        [collect(2:n); collect(1:(n-1))], 
        ones(2*(n-1)))
end

"""
    graph = star_graph(n)
"""
star_graph(n) = sparse(star_graph_ijv(n))

"""
    ijv = star_graph_ijv(n::Int64)
"""
function star_graph_ijv(n::Integer)
    IJV(n, 2*(n-1),
        [ones(Int, n-1) ; collect(2:n)], 
        [collect(2:n); ones(Int, n-1)], 
        ones(2*(n-1)))
end

"""
    graph = complete_graph(n)
"""
function complete_graph(n::Integer)
  return sparse(ones(n,n) - Matrix(I,n,n))
end 

"""
    ijv = complete_graph_ijv(n)
"""
complete_graph_ijv(n::Integer) = IJV(complete_graph(n))


"""
    graph = complete_bipartite_graph(n)
"""
function complete_bipartite_graph(n::Integer)
    return sparse([zeros(n,n) ones(n,n) ; ones(n,n) zeros(n,n)])
end 

"""
    ijv = complete_bipartite_graph_ijv(n)
"""
complete_bipartite_graph_ijv(n::Integer) = IJV(complete_bipartite_graph(n))


"""
    graph = ring_graph(n)
"""
ring_graph(n) = sparse(ring_graph_ijv(n))


"""
    ijv = ring_graph_ijv(n)
"""
function ring_graph_ijv(n::Integer)
    if n == 1
        return empty_graph_ijv(1)
    else
        return IJV(n, 2*n,
        [collect(1:(n-1)) ; collect(2:n); 1; n], 
        [collect(2:n); collect(1:(n-1)); n; 1], 
        ones(2*n))
    end
end    


"""
    graph = generalized_ring(n, gens)

A generalization of a ring graph.
The vertices are integers modulo n.
Two are connected if their difference is in gens.
For example,

```
generalized_ring(17, [1 5])
```
"""
generalized_ring(n::T, gens::Array{T}) where T <: Integer = 
    sparse(generalized_ring_ijv(n, gens))

function generalized_ring_ijv(n::T, gens::Array{T}) where T <: Integer
    gens = gens[mod.(gens, n) .> 0]
    if isempty(gens)
        return empty_graph_ijv(n)
    end

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
    return IJV(n, m, 1 .+ ai,1 .+ aj,ones(m))
end

"""
    graph = rand_gen_ring(n, k; verbose = false, ver=Vcur)

A random generalized ring graph of degree k.
Gens always contains 1, and the other k-1 edge types
are chosen from an exponential distribution
"""
rand_gen_ring(n::Integer, k::Integer; verbose=false, ver=Vcur) =
    sparse(rand_gen_ring_ijv(n,k,verbose=verbose, ver=ver))

function rand_gen_ring_ijv(n::Integer, k::Integer; verbose=false, ver=Vcur)

    # if any of n, 2n, 3n etc. is in gens we will have self loops
    gens = [0]
    while 0 in (gens .% n)
        gens = [1; 1 .+ ceil.(Integer,exp.(rand_ver(ver, k-1)*log(n-1)))]
    end

    if verbose
        println("gens: $(gens)")
    end            
    return generalized_ring_ijv(n, gens)
end



"""
    graph = hyperCube(d::Int64)

The d dimensional hypercube.  Has 2^d vertices and d*2^(d-1) edges.
"""
hypercube(d::Integer) = sparse(hypercube_ijv(d))

function hypercube_ijv(d::Integer)
    @assert d >= 0

    if d == 0
        return empty_graph_ijv(1)
    end

    ijvm = hypercube_ijv(d-1)

    ijv = disjoin(ijvm, ijvm)
    append!(ijv.i, [collect(1:2^(d-1)); 2^(d-1) .+ collect(1:2^(d-1))]  )
    append!(ijv.j, [2^(d-1) .+ collect(1:2^(d-1)); collect(1:2^(d-1))] )
    append!(ijv.v, ones(2^d))
    ijv.nnz += 2^d

    return ijv

end


"""
    graph = completeBinaryTree(n::Int64)

The complete binary tree on n vertices
"""
complete_binary_tree(n::Integer) = sparse(cbt_ijv(n))

function cbt_ijv(n::Integer)
    
    k = div(n-1,2)

    if 2*k+1 < n
        ii0 = Int[n-1]
        jj0 = Int[n]
    else
        ii0 = Int[]
        jj0 = Int[]
    end

    ii = [collect(1:k); collect(1:k); ii0]
    jj = [2*collect(1:k); 2*collect(1:k) .+ 1; jj0]

    return IJV(n, 2*(n-1),
        [ii;jj], [jj;ii], ones(2*length(ii)))

end

"""
    graph = grid2(n::Int64, m::Int64; isotropy=1)

An n-by-m grid graph.  iostropy is the weighting on edges in one direction.
"""
grid2(n::Integer, m::Integer; isotropy=1.0) = 
    sparse(grid2_ijv(n, m; isotropy=isotropy))

grid2_ijv(n::Integer, m::Integer; isotropy=1.0) =
    product_graph(isotropy*path_graph_ijv(n), path_graph_ijv(m))

grid2(n::Integer) = grid2(n,n)
grid2_ijv(n::Integer) = grid2_ijv(n,n)

"""
    graph = grid3(n1, n2, n3)
    graph = grid3(n)

An n1-by-n2-by-n3 grid graph.
"""
grid3(n1::Integer, n2::Integer, n3::Integer) = 
    sparse(grid3_ijv(n1, n2, n3))

grid3_ijv(n1::Integer, n2::Integer, n3::Integer) =
    product_graph(path_graph(n1), product_graph(path_graph(n2), path_graph(n3)))

grid3(n) = grid3(n,n,n)
grid3_ijv(n) = grid3_ijv(n,n,n)

#=
"""
    graph = wGrid2(n::Integer; weightGen::Function=rand)

An n by n grid with random weights. User can specify the weighting scheme.
"""
wgrid2(n::Integer; weightGen::Function=r(x->rand(x)) = 
    sparse(wgrid2_ijv(n, weightGen = weightGen))

function wgrid2_ijv(n::Integer; weightGen::Function=rand)
    gr2 = compress(grid2_ijv(n))

    # inefficient for backwards compatibility
    for i in 1:gr2.nnz
        gr2.v[i] = weightGen()
        if gr2.i[i] < gr2.j[i]
            gr2.v[i] = 0
        end
    end

    return compress(gr2 + gr2')

end
=#

"""
    graph = grid2coords(n::Int64, m::Int64)
    graph = grid2coords(n::Int64)

Coordinates for plotting the vertices of the n-by-m grid graph
"""
function grid2coords(n::Int64, m::Int64)
  x = kron(collect(1:n),ones(m))
  y = kron(ones(n),collect(1:m))
  return x, y
end # grid2coords

grid2coords(n) = grid2coords(n, n)


"""
    graph = rand_matching(n::Integer; ver=Vcur)

A random matching on n vertices
"""
rand_matching(n::Integer; ver=Vcur) = sparse(rand_matching_ijv(n, ver=ver))

function rand_matching_ijv(n::Integer; ver=Vcur)

  p = randperm_ver(ver,n)
  n1 = convert(Int64,floor(n/2))
  n2 = 2*n1

  ii = p[1:n1] 
  jj = p[(n1+1):n2]

  return IJV(n, n2, 
    [ii;jj], [jj; ii], ones(n2))

end 

"""
    graph = rand_regular(n, k; ver=Vcur)

A sum of k random matchings on n vertices
"""
rand_regular(n::Integer, k::Integer; ver=Vcur) = sparse(rand_regular_ijv(n, k, ver=ver))

function rand_regular_ijv(n::Integer, k::Integer; ver=Vcur)

    n1 = convert(Int64,floor(n/2))
    n2 = 2*n1

    ii = Array{Int64}(undef, n1*k)
    jj = Array{Int64}(undef, n1*k)

    ind = 0
    for i in 1:k
        p = randperm_ver(ver, n)   
        for j in 1:n1
            ind += 1
            ii[ind] = p[j]
            jj[ind] = p[n1+j]
        end
    end

    return IJV(n, k*n2,
        [ii;jj], [jj;ii], ones(k*n2))

end 


"""
    graph = grown_graph(n, k; ver=Vcur)

Create a graph on n vertices.
For each vertex, give it k edges to randomly chosen prior
vertices.
This is a variety of a preferential attachment graph.
"""
grown_graph(n::Integer, k::Integer; ver=Vcur) = sparse(grown_graph_ijv(n,k, ver=ver))


function grown_graph_ijv(n::Integer, k::Integer; ver=Vcur)
    ii = Int[]
    jj = Int[]

    for i = 1:k
        append!(ii, collect(2:n))
        append!(jj, ceil.(Integer,collect(1:n-1).*rand_ver(ver, n-1)))
    end

    return IJV(n, 2*k*(n-1),
        [ii;jj], [jj;ii], ones(2*k*(n-1)))
    
end # grownGraph

# used in grownGraphD
function randSet(n::Integer,k::Integer; ver=Vcur)
    if n == k
        return collect(1:n)
    elseif n < k
        error("n must be at least k")
    else

        s = sort(ceil.(Integer,n*rand_ver(ver, k)))
        good = (minimum(s[2:end]-s[1:(end-1)]) > 0)
        while good == false
            s = sort(ceil.(Integer,n*rand_ver(ver, k)))
            good = (minimum(s[2:end]-s[1:(end-1)]) > 0)
        end

        return s

    end
end

"""
    graph = grown_graph_d(n::Integer, k::Integer; ver=Vcur)

Like a grownGraph, but it forces the edges to all be distinct.
It starts out with a k+1 clique on the first k vertices
"""
grown_graph_d(n::Integer, k::Integer; ver=Vcur) = 
    sparse(grown_graph_d_ijv(n::Integer, k::Integer, ver=ver))

function grown_graph_d_ijv(n::Integer, k::Integer; ver=Vcur)
    @assert n > k > 1

    u = zeros(Int64, k*(n-k-1))
    v = zeros(Int64, k*(n-k-1))

    for i in (k+2):n
        nb = randSet(i-1,k; ver=ver)
        u[(i-k-2)*k .+ collect(1:k)] .= i
        v[(i-k-2)*k .+ collect(1:k)] .= nb
    end

    ijv = IJV(n, 2*length(u),
        [u;v], [v;u], ones(2*length(u)))

    clique = complete_graph_ijv(k+1)
    clique.n = n

    return ijv + clique
end 

"""
    graph = pref_attach(n::Int64, k::Int64, p::Float64; ver=Vcur)

A preferential attachment graph in which each vertex has k edges to those
that come before.  These are chosen with probability p to be from a random vertex,
and with probability 1-p to come from the endpoint of a random edge.
It begins with a k-clique on the first k+1 vertices.
"""
pref_attach(n::Integer, k::Integer, p::Float64; ver=Vcur) = sparse(pref_attach_ijv(n,k,p,ver=ver))

function pref_attach_ijv(n::Integer, k::Integer, p::Float64; ver=Vcur)
    @assert n > k
    if n == (k+1)
        return complete_graph_ijv(n)
    end

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
                if rand_ver(ver) < p
                    s[j] = rand_ver(ver,1:(i-1))
                else
                    s[j] = v[rand_ver(ver,1:(k*(i-1)))]
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

    w[1:(k*(k+1))] .= 0.5

    return IJV(n, length(w),
        [u;v], [v;u], [w;w])

end


"""
    graph = randperm(mat::AbstractMatrix)
            randperm(f::Expr)

Randomly permutes the vertex indices
"""
randperm(mat::AbstractMatrix) = randperm_ver(Vcur, mat)

function randperm_ver(::Type{V06}, mat::AbstractMatrix)
    perm = randperm_ver(V06, mat.n)
    return mat[perm,perm]
end

function randperm_ver(::Type{Vcur}, mat::AbstractMatrix)
    perm = randperm_ver(Vcur, mat.n)
    return mat[perm,perm]
end

randperm(f::Expr; ver=Vcur) = randperm(eval(f), ver=ver)

function randperm_ver(::Type{V06}, a::IJV)
    perm = invperm(randperm_ver(V06, a.n))  # invperm is for V06 complete_binary_tree
    return IJV(a.n, a.nnz,
        perm[a.i], perm[a.j], a.v)
end   

function randperm_ver(::Type{Vcur}, a::IJV)
    perm = randperm_ver(Vcur, a.n)  
    return IJV(a.n, a.nnz,
        perm[a.i], perm[a.j], a.v)
end   

function randperm_ver!(::Type{Vcur}, a::IJV)
    perm = randperm_ver(Vcur, a.n)  
    a.i = perm[a.i];
    a.j = perm[a.j];

end   

function randperm_ver!(::Type{V06}, a::IJV)
    perm = randperm_ver(V06, a.n)  
    a.i = perm[a.i];
    a.j = perm[a.j];

end   


randperm(a::IJV) = randperm_ver(Vcur, a)
  

"""
    graph = ErdosRenyi(n::Integer, m::Integer; ver=Vcur)

Generate a random graph on n vertices with m edges.
The actual number of edges will probably be smaller, as we sample
with replacement
"""
function ErdosRenyi(n::Integer, m::Integer; ver=Vcur)
    ai = rand_ver(ver, 1:n, m)
    aj = rand_ver(ver, 1:n, m)
    ind = (ai .!= aj)

    mat = sparse(ai[ind],aj[ind],1.0,n,n)
    mat = mat + mat'
    unweight!(mat)
    return mat
end

ErdosRenyi_ijv(n::Integer, m::Integer; ver=Vcur) = IJV(ErdosRenyi(n::Integer, m::Integer, ver=ver))


"""
    graph = ErdosRenyiCluster(n::Integer, k::Integer; ver=Vcur)

Generate an ER graph with average degree k,
and then return the largest component.
Will probably have fewer than n vertices.
If you want to add a tree to bring it back to n,
try ErdosRenyiClusterFix.
"""
function ErdosRenyiCluster(n::Integer, k::Integer; ver=Vcur)
    m = ceil(Integer,n*k/2)
    ai = rand_ver(ver, 1:n, m)
    aj = rand_ver(ver, 1:n, m)
    ind = (ai .!= aj)
    mat = sparse(ai[ind],aj[ind],1.0,n,n)
    mat = mat + mat'

    return biggestComp(mat)
end

ErdosRenyiCluster_ijv(n::Integer, k::Integer; ver=Vcur) = IJV(ErdosRenyiCluster(n, k, ver=ver))

"""
    graph = ErdosRenyiClusterFix(n::Integer, k::Integer; ver=Vcur)

Like an Erdos-Renyi cluster, but add back a tree so
it has n vertices
"""
function ErdosRenyiClusterFix(n::Integer, k::Integer; ver=Vcur)
    m1 = ErdosRenyiCluster(n, k)
    n2 = n - size(m1)[1]
    if (n2 > 0)
        m2 = complete_binary_tree(n2)
        return join_graphs(m1,m2,1,ver=ver)
    else
        return m1
    end
end

function ErdosRenyiClusterFix_ijv(n::Integer, k::Integer; ver=Vcur)
    m1 = ErdosRenyiCluster_ijv(n, k, ver=ver)
    n2 = n - m1.n
    if (n2 > 0)
        m2 = cbt_ijv(n2)
        join_graphs!(m1,m2,1,ver=ver)
    end

    return m1

end


"""
    graph = pure_random_graph(n::Integer; verbose=false, ver=Vcur)

Generate a random graph with n vertices from one of our natural distributions
"""
pure_random_graph(n::Integer; verbose=false, prefix="", ver=Vcur) = 
    sparse(pure_random_ijv(n; verbose=verbose, prefix=prefix, ver=ver))

function pure_random_ijv_v6(n::Integer; verbose=false, prefix="")

    gr = []
    wt = []

    push!(gr,:(path_graph_ijv($n)))
    push!(wt,1)

    push!(gr,:(ring_graph_ijv($n)))
    push!(wt,3)

    push!(gr,:(cbt_ijv($n)))
    push!(wt,3)

    push!(gr,:(grown_graph_ijv($n,2, ver=V06)))
    push!(wt,6)

    push!(gr,:(firstn(grid2_ijv(ceil(Integer,sqrt($n))), $n)))
    push!(wt,6)

    push!(gr,:(rand_regular_ijv($n,3, ver=V06)))
    push!(wt,6)

    push!(gr,:(ErdosRenyiClusterFix_ijv($n,2, ver=V06)))
    push!(wt,6)

    if n >= 4
        push!(gr,:(rand_gen_ring_ijv($n,4,ver=V06, verbose=$(verbose))))
        push!(wt,6)
    end

    i = sampleByWeight(wt, ver=V06)

    # make sure get a connected graph - will want to remove.
    if verbose
        println(prefix, gr[i])
    end
    ijv = eval(gr[i])

    its = 0
    while (~isConnected(sparse(ijv))) && (its < 100)
        i = sampleByWeight(wt, ver=V06)

        ijv = eval(gr[i])
        its += 1
    end
    if its == 100
        error("Getting a disconnected graph from $(gr[i])")
    end

    return ijv

end

"""
    a = pure_random_ijv(n::Integer; verbose=false, prefix="", ver=Vcur)

Chooses among path_graph, ring_graph, grid_graph, complete_binary_tree, rand_gen_ring, grown_graph and ErdosRenyiClusterFix.
It can produce a disconnected graph.
For code that always produces a connected graph (and is the same as with Julia v0.6, use pure_random_ijv_v6)
"""
function pure_random_ijv(n::Integer; verbose=false, prefix="", ver=Vcur)
    if ver == Vcur
        pure_random_ijv_v7(n, verbose=verbose, prefix=prefix)
    else
        pure_random_ijv_v6(n, verbose=verbose, prefix=prefix)
    end
end

function pure_random_ijv_v7(n::Integer; verbose=false, prefix="")

    ver = Vcur

    n >= 4 ? rmax = 37 : rmax = 31

    r = rmax*rand_ver(ver)

    if r <= 1
        ijv = path_graph_ijv(n)
        verbose && println("$(prefix) path_graph($(n))")

    elseif r <= 4
        ijv = ring_graph_ijv(n)
        verbose && println("$(prefix) ring_graph($(n))")

    elseif r <= 7
        ijv = cbt_ijv(n)
        verbose && println("$(prefix) complete_binary_tree($(n))")

    elseif r <= 13
        ijv = grown_graph_ijv(n,2, ver=ver)
        verbose && println("$(prefix) grown_graph($(n))")

    elseif r <= 19
        ijv = firstn(grid2_ijv(ceil(Integer,sqrt(n))), n)
        verbose && println("$(prefix) firstn_grid2($(n))")

    elseif r <= 25
        ijv = rand_regular_ijv(n,3, ver=ver)
        verbose && println("$(prefix) rand_regular_ijv($(n),3)")

    elseif r <= 31
        ijv = ErdosRenyiClusterFix_ijv(n,2, ver=ver)
        verbose && println("$(prefix) ErdosRenyiClusterFix_ijv($n,2)")

    else
        ijv = rand_gen_ring_ijv(n,4, verbose=verbose, ver=ver)   
        verbose && println("$(prefix) rand_gen_ring($(n), 4)")   
                

    end

    return ijv

end

"""
    ind = sampleByWeight(wt; ver=Vcur)

sample an index with probability proportional to its weight given here
"""
function sampleByWeight(wt; ver=Vcur)
    r = rand_ver(ver, 1)*sum(wt)
    findall(cumsum(wt) .> r)[1]
end

"""
    graph = semiwted_chimera(n::Integer; verbose=false, ver=Vcur)

A Chimera graph with some weights.  The weights just appear when graphs are combined.
For more interesting weights, use `wted_chimera`
"""
function semiwted_chimera(n::Integer; verbose=false, prefix="", ver=Vcur) 
    if ver == V06
        return sparse(semiwted_chimera_ijv_v6(n; verbose=verbose, prefix=prefix, ver=ver))
    else
        return sparse(semiwted_chimera_ijv(n; verbose=verbose, prefix=prefix, ver=ver))
    end
end

function semiwted_chimera(n::Integer, k::Integer; verbose=false, prefix="", ver=Vcur)
    if ver == V06
        seed_ver!(ver, 100*n+k)
        return sparse(semiwted_chimera_ijv_v6(n; verbose=verbose, prefix=prefix, ver=ver))
    end
        
    seed_ver!(ver, hash(n, hash(k)))

    return semiwted_chimera(n; verbose=verbose, prefix=prefix, ver=ver)

end

function semiwted_chimera_ijv(n::Integer; verbose=false, prefix="", ver=Vcur)
    if ver == V06
        return semiwted_chimera_ijv_v6(n, verbose=verbose, prefix=prefix, ver=ver)
    else
        return semiwted_chimera_ijv_v7(n, verbose=verbose, prefix=prefix, ver=ver)
    end
end       

function semiwted_chimera_ijv_v6(n::Integer; verbose=false, prefix="", ver=Vcur)

    if (n < 2)
        return empty_graph_ijv(1)
    end

    r = rand_ver(ver)^2

    if (n < 30) || (rand_ver(ver) < .2)

        gr = pure_random_ijv(n, verbose=verbose, prefix=prefix, ver=ver)

        return randperm_ver(ver, gr)
    end

    if (n < 200)
        # just join disjoint copies of graphs

        n1 = 10 + floor(Integer,(n-20)*rand_ver(ver))
        n2 = n - n1
        k = ceil(Integer,exp(rand_ver(ver)*log(min(n1,n2)/2)))

        if verbose
            println(prefix,"joinGraphs($(r)*chimera($(n1)),chimera($(n2)),$(k))")
        end

        pr = string(" ",prefix)
        gr = join_graphs(r*chimera_ijv(n1;verbose=verbose,prefix=pr, ver=ver),
          chimera_ijv(n2;verbose=verbose,prefix=pr, ver=ver),
          k, ver=ver)

        return randperm_ver(ver, gr)
    end

    # split with probability .7

    if (rand_ver(ver) < .7)
        n1 = ceil(Integer,10*exp(rand_ver(ver)*log(n/20)))

        n2 = n - n1
        k = floor(Integer,1+exp(rand_ver(ver)*log(min(n1,n2)/2)))

        if verbose
            println(prefix,"joinGraphs($(r)*chimera($(n1)),chimera($(n2)),$(k))")
        end

        pr = string(" ",prefix)
        gr = join_graphs(r*chimera_ijv(n1;verbose=verbose,prefix=pr, ver=ver),
          chimera_ijv(n2;verbose=verbose,prefix=pr, ver=ver),
          k, ver=ver)

        return randperm_ver(ver, gr)

    else
        n1 = floor(Integer,10*exp(rand_ver(ver)*log(n/100)))

        n2 = floor(Integer, n / n1)

        if (rand_ver(ver) < .5)

            if verbose
                println(prefix,"productGraph($(r)*chimera($(n1)),chimera($(n2)))")
            end
            pr = string(" ",prefix)
            gr = product_graph(r*chimera_ijv(n1;verbose=verbose,prefix=pr, ver=ver),
              chimera_ijv(n2;verbose=verbose,prefix=pr, ver=ver))
        else

            k = floor(Integer,1+exp(rand_ver(ver)*log(min(n1,n2)/10)))

            if verbose
                println(prefix, "generalizedNecklace($(r)*chimera($(n1)),chimera($(n2)),$(k))")
            end
            pr = string(" ",prefix)
            gr = generalized_necklace(r*chimera_ijv(n1;verbose=verbose,prefix=pr, ver=ver),
              chimera_ijv(n2;verbose=verbose,prefix=pr, ver=ver),
              k, ver=ver)

        end


        

        n3 = n - gr.n
        if (n3 > 0)

            if verbose
                println(prefix, "join_graphs!(gr($(gr.n)),chimera($(n3)),2), $(ver)")
            end

            pr = string(" ",prefix)
            join_graphs!(gr,chimera_ijv(n3;verbose=verbose,prefix=pr, ver=ver),
            2, ver=ver)

        end

        return randperm_ver(ver, gr)

    end
end

function semiwted_chimera_ijv_v7(n::Integer; verbose=false, prefix="", ver=Vcur)

    if (n < 2)
        return empty_graph_ijv(1)
    end

    gr = empty_graph(n)

    r = rand_ver(ver)^2

    pr = string(" ",prefix)

            
    if (n < 30) || (rand_ver(ver) < .2)

        gr = pure_random_ijv(n, verbose=verbose, prefix=prefix, ver=ver)

    elseif (n < 200)
        # just join disjoint copies of graphs

        n1 = 10 + floor(Integer,(n-20)*rand_ver(ver))
        n2 = n - n1
        k = ceil(Integer,exp(rand_ver(ver)*log(min(n1,n2)/2)))

        if verbose
            println(prefix,"joinGraphs($(r)*chimera($(n1)),chimera($(n2)),$(k))")
        end

        gr = join_graphs(r*chimera_ijv(n1;verbose=verbose,prefix=pr, ver=ver),
          chimera_ijv(n2;verbose=verbose,prefix=pr, ver=ver),
          k, ver=ver)

    elseif rand_ver(ver) < .7

            # split with probability .7

        
            n2 = ceil(Integer,10*exp(rand_ver(ver)*log(n/20)))

            if n2 < n/2
                n1 = n - n2
            else
                n1 = n2
                n2 = n - n2
            end

            k = floor(Integer,1+exp(rand_ver(ver)*log(min(n1,n2)/2)))

            if verbose
                println(prefix,"joinGraphs($(r)*chimera($(n1)),chimera($(n2)),$(k))")
            end

            gr = r*chimera_ijv(n1;verbose=verbose,prefix=pr, ver=ver)
            join_graphs!(gr,
                chimera_ijv(n2;verbose=verbose,prefix=pr, ver=ver),
                k, ver=ver)

        else
            n1 = floor(Integer,10*exp(rand_ver(ver)*log(n/100)))

            n2 = floor(Integer, n / n1)

            choice = rand_ver(ver)

            if choice < .45

                if verbose
                    println(prefix,"productGraph($(r)*chimera($(n1)),chimera($(n2)))")
                end

                gr = product_graph(r*chimera_ijv(n1;verbose=verbose,prefix=pr, ver=ver),
                chimera_ijv(n2;verbose=verbose,prefix=pr, ver=ver))

            elseif choice < 0.9

                k = floor(Integer,1+exp(rand_ver(ver)*log(min(n1,n2)/10)))

                if verbose
                    println(prefix, "generalizedNecklace($(r)*chimera($(n1)),chimera($(n2)),$(k))")
                end

                gr = generalized_necklace(r*chimera_ijv(n1;verbose=verbose,prefix=pr, ver=ver),
                chimera_ijv(n2;verbose=verbose,prefix=pr, ver=ver),
                k, ver=ver)

            else 

                n1 = div(n,2)
                k = ceil(Integer,exp(rand_ver(ver)*log(n/2)))

                if verbose
                    println(prefix, "two_lift(chimera($(n1), $(k)")
                end

                gr1 = chimera_ijv(n1;verbose=verbose, prefix=pr, ver=ver)
                gr = two_lift(gr1, k)

            end

        end
        

    n3 = n - gr.n
    if (n3 > 0)

        if verbose
            println(prefix, "join_graphs!(gr($(gr.n)),chimera($(n3)),2), $(ver)")
        end

        join_graphs!(gr,chimera_ijv(n3;verbose=verbose,prefix=pr, ver=ver),
        2, ver=ver)
    end

    randperm_ver!(ver, gr)
    return gr

end




"""
    graph = chimera(n::Integer; verbose=false, ver=Vcur)

Builds a chimeric graph on n vertices.
The components come from pureRandomGraph,
connected by joinGraphs, productGraph and generalizedNecklace
"""
function chimera(n::Integer; verbose=false, prefix="", ver=Vcur) 

    gr = sparse(chimera_ijv(n, verbose=verbose, prefix=prefix, ver=ver))

    if ver != V06
        k = 1+floor(Integer,-log.(rand_ver(ver))/2)

        if k > 1
            if verbose
                println(prefix, "thicken($(k))")    
            end

            gr = thicken(gr, k)
        end
    end

    return gr
end


function chimera_ijv(n::Integer; verbose=false, prefix="", ver=Vcur)

    ijv = semiwted_chimera_ijv(n; verbose=verbose, prefix=prefix, ver=ver)

    if ver == V06
        gr = sparse(ijv)
        unweight!(gr)
        ijv = IJV(gr)
    end

    return ijv

end


"""
    graph = chimera(n::Integer, k::Integer; verbose=false, ver=Vcur)

Builds the kth chimeric graph on n vertices.
It does this by resetting the random number generator seed.
It should captute the state of the generator before that and then
return it, but it does not yet.
"""
function chimera(n::Integer, k::Integer; verbose=false, prefix="", ver=Vcur)
    if ver == V06
        seed_ver!(ver, 100*n+k)
    else
        seed_ver!(ver, hash(n, hash(k)))
    end

    g = chimera(n; verbose=verbose, prefix=prefix, ver=ver)
    return g
end

"""
Wrapper for boundary chimera. 
It returns an SDDM matrix.
"""
function bndry_chimera(n::Integer, k::Integer; verbose=false, prefix="", ver=Vcur)
    a = chimera(n, k; verbose=verbose, prefix=prefix, ver=ver)
    L = lap(a)
    int = setdiff(1:n,1:ceil(n^(1/3)):n)
    M = L[int, int]
    return M
end

"""
Wrapper for unweighted chimera.
"""
function uni_chimera(n::Integer, k::Integer; verbose=false, prefix="", ver=Vcur)
    a = chimera(n, k; verbose=verbose, prefix=prefix, ver=ver)
    unweight!(a)
    return a 
end

"""
Wrapper for unweighted boundary chimera.
It returns an SDDM matrix.
"""
function uni_bndry_chimera(n::Integer, k::Integer; verbose=false, prefix="", ver=Vcur)
    a = chimera(n, k; verbose=verbose, prefix=prefix, ver=ver)
    unweight!(a)
    L = lap(a)
    int = setdiff(1:n,1:ceil(n^(1/3)):n)
    M = L[int, int]
    return M
end


"""
    graph = randWeight(graph; ver=Vcur)

Applies one of a number of random weighting schemes to the edges of the graph
"""
function rand_weight(a; ver=Vcur)

    if (rand_ver(ver) < .2)
        return a
    else
        return rand_weight_sub(a, ver=ver)
    end
end


function rand_weight_sub(a; ver=Vcur)

    n = a.n
    (ai,aj) = findnz(a)
    m = length(ai)

    # potentials or edge-based

    if (rand_ver(ver) < .3)
        w = rand_ver(ver, m)

    else
        v = randn_ver(ver, a.n)

        # mult by matrix ?
        if (rand_ver(ver) < .5)

            invdeg = sparse(Diagonal(1 ./(a*ones(size(a)[1]))))
            if (rand_ver(ver) < .5)
                for i in 1:10
                    v = a * (invdeg * v)
                    v = v .- mean(v)
                end
            else

                for i in 1:10
                    v = v - a * (invdeg * v)
                    v = v .- mean(v)
                end
            end
        end

        w = abs.(v[ai] - v[aj])

    end

    # reciprocate or not?

    w[w.==0] .= 1
    w[isnan.(w)] .= 1

    if (rand_ver(ver) < .5)

        w = 1 ./w
    end

    w = w / mean(w)

    ar = sparse(ai,aj,w,n,n)
    ar = ar + ar';
    return ar
end

"""
    graph = wted_chimera(n::Integer, k::Integer; verbose=false, ver=Vcur)

Builds the kth wted chimeric graph on n vertices.
It does this by resetting the random number generator seed.
It should captute the state of the generator before that and then
return it, but it does not yet.
"""
function wted_chimera(n::Integer, k::Integer; verbose=false, ver=Vcur)
    if ver == V06
        seed_ver!(ver, 100*n+k)
    else
        seed_ver!(ver, hash(n, hash(k)))
    end
    g = wted_chimera(n; verbose=verbose, ver=ver)
    return g
end


"""
    graph = wted_chimera(n::Integer)

Generate a chimera, and then apply a random weighting scheme
"""
function wted_chimera(n::Integer; verbose=false, ver=Vcur)

    gr = semiwted_chimera(n; verbose=verbose, ver=ver)


    if ver != V06
        k = 1+floor(Integer,-log.(rand_ver(ver))/2)

        if k > 1
            if verbose
                println("thicken($(k))")    
            end

            gr = thicken(gr, k)
        end
    end

    return rand_weight(gr, ver=ver)
end

function star_join(a, k)
    if !issparse(a)
        a = sparse(a)
    end

    n = size(a,1)

    anew = kron(I(k),a)
    ai, aj, av = findnz(anew)

    newv = k*n+1
    nbrs = collect(1:k)*n

    append!(ai,nbrs)
    append!(aj,newv*ones(k))
    append!(av,ones(k))

    append!(aj,nbrs)
    append!(ai,newv*ones(k))
    append!(av,ones(k))

    return sparse(ai,aj,av)
end