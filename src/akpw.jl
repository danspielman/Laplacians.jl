#=
for the unwted case first:

do a bfs from an initial vertex.
label vertices 0 initially, and eventually by component number.

as bfs, keep a boundary list.

when finish making the cluster, go down the boundary list of find the next vertex.

issue: some vertices could get added to be boundary list multiple times.
  there are two ways to deal with this:

1. start with all nodes labeled -1
when a node reaches a boundary, label it 0.
when in clusters, get positive numbers.
so, we never put in boundary twice.

we will need a reusable queue to keep track of current boundaries.

So, will be two queues: oldboundaries and currentboundary

go until bdry to vol ratio goes down enough
to keep track of vol and bdry:
  every time we add a new node to the component,
  we traverse its edges.
  each edge pointing in increases vol, and decreases bdry (was prev bdry)
  each edge pointing out increases bdry
 sounds simple enough

akpwU(graph, startNode)

init
  create boundary queue
  fastQueue has
    fixed initial size
    the queue
    endPtr
    curPtr

    cur is where are,
    endPtr is the last item

  init the comp vector to -1
  0 for bdry
  others for comp number

then, need a grow op.
  it gets to ref comp, both queues, and the graph

  it should be called with a seed node.
  So, before call, is a while loop to check on bdry.

=#


type fastQueue
    q::Array{Int64,1}
    n::Int64
    curPtr::Int64
    endPtr::Int64
end
fastQueue(n::Int) = fastQueue(zeros(Int64,n), n, 1, 0)

hasMore(fq::fastQueue) = fq.curPtr <= fq.endPtr

function push!(fq::fastQueue, i::Int64)
    @assert fq.endPtr < fq.n
    
    fq.endPtr = fq.endPtr + 1
    fq.q[fq.endPtr] = i
       
end

function pull!(fp::fastQueue)
    @assert hasMore(fq)
    
    i = fq.q[fq.curPtr]
    fq.curPtr += 1
    
    return i
end

function reset!(fq::fastQueue)
    fq.curPtr = 1
    fq.endPtr = 0
end