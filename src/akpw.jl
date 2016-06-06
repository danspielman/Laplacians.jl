#=
for the unwted case first:

do a bfs from an initial vertex.
label vertices 0 initially, and eventually by component number.


=#

#========================================================================

     DATA STRUCTURES

=#


type fastQueue
    q::Array{Int64,1}
    n::Int64
    curPtr::Int64
    endPtr::Int64
end

fastQueue(n::Int) = fastQueue(zeros(Int64,n), n, 1, 0)

hasMore(fq::fastQueue) = fq.curPtr <= fq.endPtr

import Base.push!

function push!(fq::fastQueue, i)
    @assert fq.endPtr < fq.n
    
    fq.endPtr = fq.endPtr + 1
    fq.q[fq.endPtr] = i
       
end

function pull!(fq::fastQueue)
    @assert hasMore(fq)
    
    i = fq.q[fq.curPtr]
    fq.curPtr += 1
    
    return i
end

function reset!(fq::fastQueue)
    fq.curPtr = 1
    fq.endPtr = 0
end

type fastPairQueue
    q::Array{Int64,2}
    n::Int64
    curPtr::Int64
    endPtr::Int64
end

fastPairQueue(n::Int) = fastPairQueue(zeros(Int64,n,2), n, 1, 0)

hasMore(fq::fastPairQueue) = fq.curPtr <= fq.endPtr

function push!(fq::fastPairQueue, i, j)
    @assert fq.endPtr < fq.n
    
    fq.endPtr = fq.endPtr + 1
    fq.q[fq.endPtr,1] = i
    fq.q[fq.endPtr,2] = j
       
end

function pull!(fq::fastPairQueue)
    @assert hasMore(fq)
    
    i1 = fq.q[fq.curPtr,1]
    i2 = fq.q[fq.curPtr,2]
    fq.curPtr += 1
    
    return (i1,i2)
end

function reset!(fq::fastPairQueue)
    fq.curPtr = 1
    fq.endPtr = 0
end




#========================================================================

              MAIN FUNC HERE

=#


function akpwU(graph, xfac::Float64)
    akpwU(graph, n->xfac)
end


function akpwU(graph, xfac::Function)
    n = size(graph,1)

    tre = akpwUsub(graph, xfac)

    (ai,aj,av) = findnz(graph)
    tree = sparse([ai[tre];aj[tre]], [aj[tre];ai[tre]], [av[tre];av[tre]])

    return tree

end

akpwU(graph) = akpwU(graph, x->(1/(2*log(x))))


function akpwUsub(graph, xfac::Function)
    n = size(graph,1)
    m = nnz(graph)
    
    comp = zeros(Int, n)
    ncomps = 0
    
    # the indices into (ai,aj) = findnz(a) of the edges in the tree
    treeEdges = fastQueue(n-1)
    
    curBdry = fastQueue(n)
    
    # start from a random vertex
    seedlist = randperm(n)
    
    # allocate reusable queue, but only used inside subroutines
    thisBdry = fastPairQueue(m)

    for seed in seedlist

        if comp[seed] <= 0
            ncomps += 1
            bfsFromSeed(graph, seed, ncomps, comp, treeEdges, thisBdry, xfac(n)) 
            reset!(thisBdry)
        end
        
    end
     
    #println(comp)
    
    #=  THIS CODE CHECKS THAT THE FOREST SPANS THE COMPONENTS 
    @assert 2*n - nnz(tree) - 2*maximum(comp) == 0
    for i in 1:maximum(comp)
        ind = find(comp .== i)
        tri = tree[ind,ind]
        @assert isConnected(tri)
    end
    =#

    tre = treeEdges.q[1:treeEdges.endPtr]


    if maximum(comp) > 1

        cGraph = compGraphU(graph, comp)

        edgeMap = cGraph.nzval
        cGraph.nzval = ones(length(edgeMap))

        if (nnz(cGraph) > 0)
            ctre = akpwUsub(cGraph, xfac)

            sube = edgeMap[ctre]

            tre = [tre;sube]
        end
        
    end
    
        
    return tre
    
end




function bfsFromSeed(graph, seed::Int, ncomps::Int, comp, 
    treeEdges::fastQueue, thisBdry::fastPairQueue, xfac::Float64)
   
    bdry = 0
    vol = 0
    
    comp[seed] = ncomps
    for ind in graph.colptr[seed]:(graph.colptr[seed+1]-1)
        nbr = graph.rowval[ind]
        if comp[nbr] == 0
            
            push!(thisBdry,ind,nbr)
            bdry += 1
            vol += 1 # if using sum of degrees CHECK ON THIS

        end
    end

    while (bdry > xfac*vol) && hasMore(thisBdry)
        (edge,node) = pull!(thisBdry)

        #println([edge node bdry vol])
        
        # if (comp[node] != ncomps)
        if (comp[node] == 0)

            comp[node] = ncomps 
            push!(treeEdges,edge)

            for ind in graph.colptr[node]:(graph.colptr[node+1]-1)
                nbr = graph.rowval[ind]
                if comp[nbr] == ncomps
                    bdry -= 1
                    vol += 1

                elseif comp[nbr] == 0

                    push!(thisBdry,ind,nbr)  # issue: nodes pop up many times
                    bdry += 1
                    vol += 1 # if using sum of degrees CHECK ON THIS

                end
            end
        end
    end

end


# unweighted version      
function compGraphU(graph, comp)

    m = nnz(graph)
    (ai,aj) = findnz(graph)

    aic = comp[ai]
    ajc = comp[aj]

    nc = maximum(comp)

    combine(a,b) = a

    aind = collect(1:m)

    # remove self-loops
    for i in 1:m
        if aic[i] == ajc[i]
            aind[i] = 0
        end
    end

    cGraph = sparse(aic,ajc,aind,nc,nc,combine)

    return cGraph

end


function akpwish(graph, xfac::Float64)
    akpwish(graph, n->xfac)
end


function akpwish(graph, xfac::Function)
    n = size(graph,1)

    tre = akpwSub(graph, xfac)

    (ai,aj,av) = findnz(graph)
    tree = sparse([ai[tre];aj[tre]], [aj[tre];ai[tre]], [av[tre];av[tre]])

    return tree

end

akpwish(graph) = akpwish(graph, x->(1/(2*log(x))))



function akpwSub(graph, xfac::Function)
    n = size(graph,1)
    m = nnz(graph)

    (_,aj,_) = findnz(graph)
    
    comp = zeros(Int, n)
    ncomps = 0
    
    # the indices into (ai,aj) = findnz(a) of the edges in the tree
    treeEdges = fastQueue(n-1)
    
    curBdry = fastQueue(n)
    
    # start from a random vertex
    # seedlist = randperm(n)

    # order the edges, keeping only those within xfac of largest
    bigEdges = find(graph.nzval .> xfac(n)*maximum(graph.nzval))
    bigOrder = sortperm(graph.nzval[bigEdges],rev=true)
    
    edgeOrder = bigEdges[bigOrder]
    
    
    # allocate reusable queue, but only used inside subroutines
    thisBdry = fastPairQueue(m)

    dists = Inf*ones(Float64, n)

    # for seed in seedlist
    for edge in edgeOrder
        edgeu = aj[edge]
        edgev = graph.rowval[edge]
        if (comp[edgeu] == 0) && (comp[edgev] == 0)

            seed = edgeu
            
            ncomps += 1
            dijkstraFromSeed(graph, aj, seed, ncomps, comp, dists, treeEdges, xfac(n)) 
            reset!(thisBdry)
        end
        
    end


    # clean up: make singletons their own comps...
    for i in 1:n
        if comp[i] == 0
            ncomps += 1
            comp[i] = ncomps
        end
    end
    

    
    #println(comp)
    
    tre = treeEdges.q[1:treeEdges.endPtr]
     
    #=  THIS CODE CHECKS THAT THE FOREST SPANS THE COMPONENTS 
    AND, computes the stretch of each 

    (ai,aj,av) = findnz(graph)
    tree = sparse([ai[tre];aj[tre]], [aj[tre];ai[tre]], [av[tre];av[tre]],n,n)

    @assert 2*n - nnz(tree) - 2*maximum(comp) == 0
    for i in 1:maximum(comp)
        ind = find(comp .== i)
        # println(ind)
        tri = tree[ind,ind]
        @assert isConnected(tri)
    end

    for i in 1:maximum(comp)
        ind = find(comp .== i)
        tri = tree[ind,ind]
        gri = a[ind,ind]
        println([i sum(compStretches(tri,gri))])
        println(ind)
    end

    =#

    if maximum(comp) > 1

        #println(comp)
        
        cGraph, edgeMap = compGraph(graph, comp)

        #println(cGraph)
        
        #=
        cGraph = compGraphU(graph, comp)
        edgeMap = cGraph.nzval
        cGraph.nzval = ones(length(edgeMap))
        =#
        
        if (nnz(cGraph) > 0)
            ctre = akpwSub(cGraph, xfac)

            sube = edgeMap[ctre]

            tre = [tre;sube]
        end
        
    end
        
    return tre
    
end




# grow shortest path tree from the seed
# will need to make the heap reusable
function dijkstraFromSeed(graph, aj::Array{Int64,1}, seed::Int, ncomps::Int, comp, 
    dists::Array{Float64,1}, treeEdges::fastQueue, xfac::Float64)

    ai = graph.rowval
    av = graph.nzval
    
    bdry = 0
    vol = 0

    heap = Collections.PriorityQueue(Int64, Float64)

    dists[seed] = 0
    comp[seed] = ncomps


   #=
   issue: how to keep track of the edges used.
    shortestPaths does it by keeping a parent pointer for every node.
    randishPrim does it by popping.
   =#

    for ind in graph.colptr[seed]:(graph.colptr[seed+1]-1)
        nbr = graph.rowval[ind]
        if comp[nbr] == 0

            wt = graph.nzval[ind]
            Collections.enqueue!(heap, ind, 1/wt)
            bdry += wt
            vol += wt
        end
        
    end

    while (bdry > xfac*vol) && (length(heap) > 0)

        edge = Collections.dequeue!(heap)
        node = aj[edge]
        from = ai[edge]
        if comp[node] == ncomps
            node = ai[edge]
            from = aj[edge]
        end


        # println([edge node bdry vol length(heap)])
        
        if (comp[node] == 0)

            comp[node] = ncomps
            # println([edge ai[edge] aj[edge] ncomps])
            push!(treeEdges,edge)
            dist = dists[from] + 1/graph.nzval[edge]
            dists[node] = dist

            #println(dist)

            for ind in graph.colptr[node]:(graph.colptr[node+1]-1)
                nbr = graph.rowval[ind]

                if comp[nbr] == ncomps
                    wt = graph.nzval[ind]

                    bdry -= wt
                    vol += wt

                elseif comp[nbr] == 0

                    wt = graph.nzval[ind]
                    newdist = dist + 1/wt
                    
                    bdry += wt
                    vol += wt

                    # not clear why the effect of this line is different from what follows..except for having a smaller heap.
                    # so, there is a bug somewhere!
                    Collections.enqueue!(heap,ind,newdist)
                    #=
                    if newdist < dists[nbr]
                        dists[nbr] = newdist
                        Collections.enqueue!(heap,ind,newdist)
                    end
                    =#
                end
            end
        end
    end

end


# this combines all the vertices in a comp together, keeping the max of their weights
# it also produces a map from edges in the smaller graph back to their reps in the original
# possible speed up: combine the third pass with the second
# possible improvement: keep track of sum of multiedges in some way
function compGraph(graph, comp)

    m = nnz(graph)
    (ai,aj,av) = findnz(graph)

    aic = comp[ai]
    ajc = comp[aj]

    nc = maximum(comp)

    aind = collect(1:m)

    # remove self-loops
    for i in 1:m
        if aic[i] == ajc[i]
            aind[i] = 0
        end
    end

    ind = find(aind)
    aic = aic[ind]
    ajc = ajc[ind]
    av = av[ind]
    aind = aind[ind]

    # println([aic ajc av aind])
    
    cGraph, edgeMap = combineMultiG(aic,ajc,av,aind)

    return cGraph, edgeMap

end

# question is how to combine.  right now use max of wts, but sum might be reasonable too
# it carries around with it aind--which is a map on edges
function combineMultiG{Ti,Tv}(ai::Array{Ti,1}, aj::Array{Ti,1}, av::Array{Tv,1},  aind::Array{Ti,1})

    numnz = length(ai)

    n = maximum([maximum(ai) maximum(aj)]) 
    
    deg = zeros(Ti, n) 

    ptr = 1
    for k in 1:numnz
        deg[ai[k]] += 1
    end

    I1 = Array(Ti, numnz)
    J1 = Array(Ti, numnz)
    aind1 = Array(Ti, numnz)
    V1 = zeros(Tv, numnz)

    cumdeg = cumsum(deg)
    col = [1;cumdeg+1]

    cumdeg1 = copy(cumdeg)
    
    for i in numnz:-1:1
        ptr = cumdeg1[ai[i]]
        cumdeg1[ai[i]] -= 1
        I1[ptr] = ai[i]
        J1[ptr] = aj[i]
        V1[ptr] = av[i]
        aind1[ptr] = aind[i]
    end


    
    I2 = Array(Ti, numnz)
    J2 = Array(Ti, numnz)
    aind2 = Array(Ti, numnz)
    V2 = Array(Tv, numnz)

    for i in numnz:-1:1
        ptr = cumdeg[J1[i]]
        cumdeg[J1[i]] -= 1
        I2[ptr] = I1[i]
        J2[ptr] = J1[i]
        aind2[ptr] = aind1[i]
        V2[ptr] = V1[i]
    end

    # now, the list is sorted by J, and within that by i
    # so compress it
    
    I3 = Array(Ti, numnz)
    J3 = Array(Ti, numnz)
    aind3 = Array(Ti, numnz)
    V3 = Array(Tv, numnz)
    
    ptr = 0
    i2old = I2[1]
    j2old = J2[1]
    v2 = 0
    ind3 = 0

    for i in 1:numnz
        i2 = I2[i]
        j2 = J2[i]
        
        
        if  (i2 == i2old) && (j2 == j2old)
            if V2[i] > v2
                v2 = V2[i]
                ind3 = aind2[i]
            end
        else 
            
            
            ptr += 1
            I3[ptr] = i2old
            J3[ptr] = j2old
            aind3[ptr] = ind3
            V3[ptr] = v2
            i2old = i2
            j2old = j2
            v2 = V2[i]
            ind3 = aind2[i]
        end

    end
    
    ptr += 1
    I3[ptr] = i2old
    J3[ptr] = j2old
    aind3[ptr] = ind3
    V3[ptr] = v2

    J3 = J3[1:ptr]
    I3 = I3[1:ptr]
    V3 = V3[1:ptr]
    aind3 = aind3[1:ptr]
    
    return sparse(J3, I3, V3, n , n), aind3 

end


