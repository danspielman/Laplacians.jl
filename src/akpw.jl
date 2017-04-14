#=
akpw.jl, part of Laplacians.jl.
By Daniel A. Spielman
Inspired by the paper
A Graph-Theoretic Game and Its Application to the k-Server Problem,
SIAM J. Comput., 24(1), 78–100,
by Noga Alon, Richard M. Karp, David Peleg, and Douglas West

This code is somewhat more aggressive than the algorithm suggested in that paper.

=#

using DataStructures

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


type reusableIntMap
    q::fastQueue
    map::Array{Int,1}
end

reusableIntMap(n::Int) = reusableIntMap(fastQueue(n), zeros(n))

# note: does not allow for resetting
function set!(rim::reusableIntMap, from::Int, to::Int)
    if rim.map[from] == 0
        push!(rim.q, from)
        rim.map[from] = to
        return true
    else
        return false
    end
end

import Base.getindex
getindex(rim::reusableIntMap, from::Int) = rim.map[from]

function reset!(rim::reusableIntMap)
    while hasMore(rim.q) > 0
        from = pull!(rim.q)
        rim.map[from] = 0
    end
    reset!(rim.q)
end





immutable IJVind
    i::Int64
    j::Int64
    v::Float64
    ind::Int64
end

import Base.isless
isless(x::IJVind, y::IJVind) = x.v < y.v

# requires sorting on i and j, with j primary
type IJVindGraph
    list::Array{IJVind,1}
    colptr::Array{Int64,1}
end

import Base.getindex
getindex(G::IJVindGraph, i::Int) = G.list[i]

function printijv(ijv::IJVind)
    println(ijv.i, " ", ijv.j, " ", ijv.v, " ", ijv.ind )
end
function printijv(list::Array{IJVind,1})
    for ijv in list
        printijv(ijv)
    end
end




#========================================================================

              Unweighted Algorithm Here

=#


"""
    tree = akpwU(graph)

Computes a low stretch spanning tree of an unweighted `graph`, and returns it as a graph.
"""
function akpwU(graph)
    n = size(graph,1)

    tre = akpwUsub(graph)

    (ai,aj,av) = findnz(graph)
    tree = sparse([ai[tre];aj[tre]], [aj[tre];ai[tre]], [av[tre];av[tre]])

    return tree

end


function akpwUsub(graph)
    n = size(graph,1)
    m = nnz(graph)
    
    comp = zeros(Int, n)
    ncomps = 0
    
    # the indices into (ai,aj) = findnz(a) of the edges in the tree
    treeEdges = fastQueue(n-1)
    
    # start from a random vertex
    seedlist = randperm(n)
    
    # allocate reusable queue, but only used inside subroutines
    thisBdry = fastPairQueue(m)

    xf = 1/(2*log(n))

    for seed in seedlist

        if comp[seed] <= 0
            ncomps += 1
            bfsFromSeed(graph, seed, ncomps, comp, treeEdges, thisBdry, xf) 
            reset!(thisBdry)
        end
        
    end
     
    tre = treeEdges.q[1:treeEdges.endPtr]

    if maximum(comp) > 1

        cGraph = compGraphU(graph, comp)

        edgeMap = copy(cGraph.nzval)
        for i in 1:length(edgeMap)
            cGraph.nzval[i] = 1
        end
        
        if (nnz(cGraph) > 0)
            ctre = akpwUsub(cGraph)

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
            vol += 1 

        end
    end

    while (bdry > xfac*vol) && hasMore(thisBdry)
        (edge,node) = pull!(thisBdry)

        if (comp[node] == 0)

            comp[node] = ncomps 
            push!(treeEdges,edge)

            for ind in graph.colptr[node]:(graph.colptr[node+1]-1)
                nbr = graph.rowval[ind]
                if comp[nbr] == ncomps
                    bdry -= 1
                    vol += 1

                elseif comp[nbr] == 0

                    push!(thisBdry,ind,nbr) 
                    bdry += 1
                    vol += 1 # if using sum of degrees CHECK ON THIS

                end
            end
        end
    end

end


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


#========================================================================

              Weighted Algorithm 

=#


"""
    tree = akpw(graph; ver=0)

Computes a low stretch spanning tree of `graph`, and returns it as a graph.
The default version is 0.  In event of emergency, one can try `ver=2`.  It is usually slower, but might have slightly better stretch.
"""
function akpw(graph; ver=0)
    n = size(graph,1)

    if ver == 0
        tre = akpwSub5(graph)
    elseif ver == 2
        tre = akpwSub2(graph)
    end

    (ai,aj,av) = findnz(graph)
    tree = sparse([ai[tre];aj[tre]], [aj[tre];ai[tre]], [av[tre];av[tre]])

    return tree

end


function akpwSub5(graph)
    n = size(graph,1)

    (ai,aj,av) = findnz(graph)

    m = length(ai)

    origList = sort(IJVindList(ai,aj,av),rev=true)

    # the indices into (ai,aj) = findnz(a) of the edges in the tree
    treeEdges = fastQueue(n-1)

    nameMap = IntDisjointSets(n)

    xf = 1/(2*log(n))

    # figure out how far to go down the list : should be edges 1:last
    maxv = origList[1].v
    last = 2
    targ::Float64 = xf*maxv
    while (last <= m) && (origList[last].v > targ) 
        last += 1
    end
    last -= 1

    curIJVind = origList[1:last]

    rim = reusableIntMap(n)

    nverts = compressIndices!(curIJVind::Array{IJVind,1}, rim::reusableIntMap)

    # this could be a lot of the time used
    curIJVind = sortIJVind(curIJVind)
    # curIJVind = compress(curIJVind)

    while (last <= m) && (nverts > 1) 

        prevTreePtr = treeEdges.endPtr

        nleft = n - treeEdges.endPtr
        xf = 1/(2*log(nleft))

        cluster!(curIJVind, treeEdges, xf) 


        for i in (prevTreePtr+1):treeEdges.endPtr
            edgeind = treeEdges.q[i]
            ainame = find_root(nameMap, ai[edgeind])
            ajname = find_root(nameMap, aj[edgeind])
            if (ainame < ajname)
                union!(nameMap, ainame,ajname)
            else
                union!(nameMap, ajname,ainame)
            end
        end

        # make the new curList, by applying NameMap to cur ind
        # and, find the max wt edge between clusters
        # remove self loops as go

        newIJVind = Array(IJVind,0)
        maxv = 0
        for i in 1:length(curIJVind)
            ijv = curIJVind[i]
            ind = ijv.ind
            namei = find_root(nameMap,ai[ind])
            namej = find_root(nameMap,aj[ind])
            if (namei != namej)
                if (ijv.v > maxv)
                    maxv = ijv.v
                end
                push!(newIJVind, IJVind(namei, namej, ijv.v, ijv.ind))
            end
        end

        last += 1

        while (last <= m) && (origList[last].v > xf*maxv)
            ijv = origList[last]
            
            namei = find_root(nameMap,ijv.i)
            namej = find_root(nameMap,ijv.j)
            if namei != namej
                if maxv == 0
                    maxv = ijv.v
                end
                push!(newIJVind, IJVind(namei, namej,  ijv.v, ijv.ind))
            end
            last += 1
        end
        last -= 1

        # would it be better to do these in bulk?
        # append!(curList,origList[prevlast:last])

                        
        nverts = compressIndices!(newIJVind::Array{IJVind,1}, rim::reusableIntMap)

        if nverts > 1
            # this can also be a big bunch of time
            curIJVind = sortIJVind(newIJVind)
            curIJVind = compress(curIJVind)
        end
        
        
    end

    tre = treeEdges.q[1:treeEdges.endPtr]
    return tre
    
end




# based on counting sort: is stable.  exploits symmetry, can produce multiedges
function IJVindGraph(inList::Array{IJVind,1})

    Ti = Int64
    Tv = Float64
    
    numnz = length(inList)
    n = 0
    for i in 1:numnz
        if inList[i].i > n
            n = inList[i].i
        end
    end
    
    deg = zeros(Ti, n) 

    ptr = 1
    for i in 1:numnz
        deg[inList[i].i] += 1
    end


    list1 = Array(IJVind, numnz)
    
    cumdeg = cumsum(deg)
    colptr = [1;cumdeg+1]

    cumdeg1 = copy(cumdeg)

    for i in numnz:-1:1
        thisi = inList[i].i
        ptr = cumdeg1[thisi]
        cumdeg1[thisi] -= 1
        list1[ptr] = inList[i]
    end

    cumdeg1 = copy(cumdeg)

    list2 = Array(IJVind, numnz)

    for i in numnz:-1:1
        thisj = list1[i].j
        ptr = cumdeg1[thisj]
        cumdeg1[thisj] -= 1
        list2[ptr] = list1[i]
    end

    return IJVindGraph(list2, colptr)

end


# based on counting sort: is stable.  exploits symmetry, can produce multiedges
# tried to improve this in sortIJVind3, but it was not faster
function sortIJVind(inList::Array{IJVind,1})

    Ti = Int64
    Tv = Float64
    
    numnz = length(inList)
    n = 0
    for i in 1:numnz
        if inList[i].i > n
            n = inList[i].i
        end
    end
    
    deg = zeros(Ti, n) 

    ptr = 1
    for i in 1:numnz
        deg[inList[i].i] += 1
    end


    list1 = Array(IJVind, numnz)
    
    cumdeg = cumsum(deg)
    cumdeg1 = copy(cumdeg)

    for i in numnz:-1:1
        thisi = inList[i].i
        ptr = cumdeg1[thisi]
        cumdeg1[thisi] -= 1
        list1[ptr] = inList[i]
    end

    cumdeg1 = copy(cumdeg)

    list2 = Array(IJVind, numnz)

    for i in numnz:-1:1
        thisj = list1[i].j
        ptr = cumdeg1[thisj]
        cumdeg1[thisj] -= 1
        list2[ptr] = list1[i]
    end

    return list2

end


# combine multiedges by keeping the one of max weight.
# must be sorted for this to work
function compress(inList::Array{IJVind,1})

    outlist = Array(IJVind,0)

    ijv = inList[1]

    iold = ijv.i
    jold = ijv.j
    v = ijv.v
    ind = ijv.ind

    for i in 2:length(inList)
        ijv = inList[i]

        if  (ijv.i == iold) && (ijv.j == jold)
            if ijv.v > v
                v = ijv.v
                ind = ijv.ind
            end
        else 

            push!(outlist, IJVind(iold, jold, v, ind))
            iold = ijv.i
            jold = ijv.j
            v = ijv.v
            ind = ijv.ind
        end
    end

    push!(outlist, IJVind(iold, jold, v, ind))

    return outlist

end


# convert (ai,aj,av) to an array of IJVind entries
function IJVindList{Tv,Ti}(ai::Array{Ti,1},aj::Array{Ti,1},av::Array{Tv,1})

    m = length(ai)
    origList = Array(IJVind,m)
    for i in 1:m
        origList[i] = IJVind(ai[i],aj[i],av[i],i)
    end

    return origList

end

# this maps the indices to consecutive integers starting at 1
function compressIndices!(curIJVind::Array{IJVind,1}, rim::reusableIntMap)
    nverts = 1

    #println("comp in ")
    #printijv(curIJVind)

    for ijv in curIJVind
        if set!(rim, ijv.j, nverts)
            nverts += 1
        end
    end

    for i in 1:length(curIJVind)
        ijv = curIJVind[i]
        curIJVind[i] = IJVind(rim[ijv.i], rim[ijv.j], ijv.v, ijv.ind)
    end

    reset!(rim)

    return (nverts-1)
end

        

    

function cluster!(curIJVind, treeEdges, xf) 

    n = curIJVind[end].j
    @assert n >= curIJVind[end].i

    # create colptr
    deg = zeros(Int, n) 

    ptr = 1
    for ijv in curIJVind
        deg[ijv.j] += 1
    end
    cumdeg = cumsum(deg)
    colptr = [1;cumdeg+1]


    ijvGraph = IJVindGraph(curIJVind, colptr) 
    
    comp = zeros(Int,n)

    ncomps = 0


    for seed in 1:n
        ijvind = curIJVind[seed]
        edgeu = ijvind.i
        edgev = ijvind.j

        #println("seed ", ijvind.v)
        
        if (comp[edgeu] == 0) && (comp[edgev] == 0)

            seed = edgeu
            ncomps += 1
            dijkstraFromSeed(ijvGraph, seed, ncomps, comp, treeEdges, xf)
            
        end
        
    end

end


    


immutable HeapEntry
    node::Int64
    edge::Int64
    dist::Float64
end

import Base.isless
isless(x::HeapEntry, y::HeapEntry) = x.dist < y.dist


# grow shortest path tree from the seed
# might want to make the heap reusable
# might want to store both edge and vertex on the heap, too
function dijkstraFromSeed(ijvGraph::IJVindGraph, seed::Int, ncomps::Int, comp, 
    treeEdges::fastQueue, xfac::Float64)

    bdry = 0
    vol = 0

    heap = Array(HeapEntry, 0)

    comp[seed] = ncomps

    for ind in ijvGraph.colptr[seed]:(ijvGraph.colptr[seed+1]-1)
        nbr = ijvGraph[ind].i 
        if comp[nbr] == 0

            wt = ijvGraph[ind].v
            Collections.heappush!(heap, HeapEntry(nbr, ijvGraph[ind].ind, 1/wt))
            bdry += wt
            vol += wt
        end
        
    end

    while (bdry > xfac*vol) && (length(heap) > 0)

        he = Collections.heappop!(heap)

        node = he.node
        
        if (comp[node] == 0)

            comp[node] = ncomps

            push!(treeEdges,he.edge)
            dist = he.dist

            for ind in ijvGraph.colptr[node]:(ijvGraph.colptr[node+1]-1)

                nbr = ijvGraph[ind].i 

                if comp[nbr] == ncomps

                    wt = ijvGraph[ind].v

                    bdry -= wt
                    vol += wt

                elseif comp[nbr] == 0

                    wt = ijvGraph[ind].v
                    newdist = dist + 1/wt
                    
                    bdry += wt
                    vol += wt

                    Collections.heappush!(heap, HeapEntry(nbr, ijvGraph[ind].ind, newdist))
                end
            end
        end
    end

end



#==============================================================

   alt variant of the code: this is a recursive version
   more closely related to the undirected one.
   
   it is a little bit slower, so we don't us it by default    
   but, its stretch is often a little lower

=#

function akpwSub2(graph)
    n = size(graph,1)


    #println("n ", n)

    (ai,aj,av) = findnz(graph)
    m = length(ai)
    
    origList = Array(IJVind,m)
    for i in 1:m
        origList[i] = IJVind(ai[i],aj[i],av[i],i)
    end
    

    origList = sort(origList,rev=true)

    comp = zeros(Int, n)
    ncomps = 0
    
    # the indices into (ai,aj) = findnz(a) of the edges in the tree
    treeEdges = fastQueue(n-1)
    
    # figure out how far to go down the list : should be edges 1:last
    maxv = origList[1].v
    last = 2
    xf = 1/(2*log(n))

    while (last <= m) && (origList[last].v > xf*maxv) 
        last += 1
    end
    last -= 1

    ijvGraph = IJVindGraph(origList[1:last]) 
    
    # for seed in seedlist
    for ijvind in origList[1:last] # ITER ON IND INSTEAD?
        edgeu = ijvind.i
        edgev = ijvind.j
        
        if (comp[edgeu] == 0) && (comp[edgev] == 0)

            seed = edgeu
            
            ncomps += 1
            dijkstraFromSeed(ijvGraph, seed, ncomps, comp, treeEdges, xf)
            
        end
        
    end



    # clean up: make singletons their own comps...
    for i in 1:n
        if comp[i] == 0
            ncomps += 1
            comp[i] = ncomps
        end
    end
    
    
    tre = treeEdges.q[1:treeEdges.endPtr]

    if maximum(comp) > 1

        
        cGraph, edgeMap = compGraph(graph, comp)

        if (nnz(cGraph) > 0)
            ctre = akpwSub2(cGraph)

            sube = edgeMap[ctre]

            tre = [tre;sube]
        end
        
    end
        
    return tre
    
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



#==============================================

   Unused Code



using Base.Order
# This was an attempt to improve on sortIJVind by sorting in place
# it does not seem to be faster in general
function sortIJVind3(inList::Array{IJVind,1})

    Ti = Int64
    Tv = Float64
    
    numnz = length(inList)
    n = 0
    for i in 1:numnz
        if inList[i].i > n
            n = inList[i].i
        end
    end
    
    deg = zeros(Ti, n) 

    ptr = 1
    for i in 1:numnz
        deg[inList[i].i] += 1
    end


    list1 = Array(IJVind, numnz)

    maxdeg = maximum(deg)
    
    cumdeg = cumsum(deg)
    colptr = [1;cumdeg+1]

    cumdeg1 = copy(cumdeg)

    for i in numnz:-1:1
        thisj = inList[i].j
        ptr = cumdeg1[thisj]
        cumdeg1[thisj] -= 1
        list1[ptr] = inList[i]
    end

    byi(x::IJVind) = x.i
    o = ord(isless, byi, false, Forward)

    for j in 1:n
        sort!(list1, colptr[j], (colptr[j+1]-1),  PartialQuickSort(colptr[j]:(colptr[j+1]-1)), o)
    end

    return list1

end
=#
