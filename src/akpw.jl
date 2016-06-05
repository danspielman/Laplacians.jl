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
            growFromSeed(graph, seed, ncomps, comp, treeEdges, thisBdry, xfac(n)) 
            reset!(thisBdry)
        end
        
    end
     
    
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

        cGraph = compGraph(graph, comp)

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




function growFromSeed(graph, seed::Int, ncomps::Int, comp, 
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
        
        if (comp[node] != ncomps)

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


      
function compGraph(graph, comp)

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

