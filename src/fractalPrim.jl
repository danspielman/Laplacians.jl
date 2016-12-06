using DataStructures
using Laplacians


immutable ValWtPair
    val::Int64
    wt::Float64
    st_pt::Int64
end

immutable StartPt
    val::Int64
    offset::Float64 # Flat amount added
end

import Base.isless

isless(e1::ValWtPair, e2::ValWtPair) = e1.wt < e2.wt

# Runs shortest path starting from each of the starting points. Returns a list
# of the new starting points. There is one new point for each original starting
# point, where for each starting point v, there is a new point u where u is the
# furthest point among the set of points closets to v in the graph.
function getNewStartPoints{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind},
                                      start_pts::Array{StartPt,1})
    n = mat.n
    (ai, aj, av) = findnz(mat)

    # Connected components of graph
    comps = IntDisjointSets(n)
        
    s = Vector{ValWtPair}()
    # Add edges
    # st_pt will hold the index of the start point.
    for i in 1:length(start_pts)
        st = start_pts[i]
        v = st.val
        for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
            wt = mat.nzval[ind]
            Base.Collections.heappush!(s, ValWtPair(ind, 1/wt, i))
        end
    end

    # Need to record each new point and the distance to that point
    f_c = Vector{StartPt}(length(start_pts))
    for i in 1:length(start_pts)
        f_c[i] = StartPt(0, 0)
    end
    while comps.ngroups != 1

      valwt = Base.Collections.heappop!(s) #hog 
      edge = valwt.val

      v = ai[edge]

      if !DataStructures.in_same_set(comps, ai[edge], aj[edge])
          DataStructures.union!(comps, ai[edge], aj[edge])
          st_pt = valwt.st_pt
          # We can assign the furthest closets of this point to be v
          # Because v will only be added if it is at least as far
          # as the previous furthest
          f_c[st_pt] = StartPt(v, valwt.wt)
          for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
              u = mat.rowval[ind]
              wt = mat.nzval[ind]
              dist = valwt.wt + 1/wt
              if !DataStructures.in_same_set(comps, v, u)
                  Base.Collections.heappush!(s, ValWtPair(ind, dist, valwt.st_pt))
              end
          end
      end
    end
    return f_c
end # getNewStartPoints

function makeStartPoints{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind},
                                    num_rounds::Int64)
    n = mat.n
     # random starting point
    v = rand(1:n)

    # Making starting points
    start = Vector{StartPt}()
    push!(start, StartPt(v, 0)) 

    # Generating starting points
    for i in 1:num_rounds
        # Furtherest closest
        f_c = getNewStartPoints(mat, start)
        # Displacement
        # OLD DISPLACEMENT CODE
        # disp = 0
        for j in 1:length(f_c)
            pt = f_c[j].val
            dist = f_c[j].offset
            adjusted_dist = sqrt(dist)
            # If the new starting point was not set, then that means that the
            # previous starting point had no points closer to it than any
            # other starting point. Thus it is a bad starting point so we
            # remove it from our list
            # OLD DISPLACEMENT CODE
            # UPDATE: removing this point from the list seems like it might not be advisable
            # REASON: in the next 'round' the some point will almost certainly be found again
            # leading to it being readded over and over
            # Code for 'displacement' commented out
            if pt == 0
                # splice!(start, j-disp)
                # disp += 1
                continue
            end
            # What choice of offset should we have?
            # Currently is distance from other start + their offset
            # OLD DISPLACEMENT CODE
            # offset = dist + start[j-disp].offset
            offset = adjusted_dist + start[j].offset
            # Currently setting weight to be round number, not sure if should
            # keep constant for starting points?
            wt = i
            push!(start, StartPt(pt, offset)) 
        end
    end
    #for st in start
    #    @printf("%s",typeof(st))
    #end

    return start
end # makeStartPoints

function makeUniformExpRVStartPoints{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind}, scale::Tval, seed::Int64 = 1)
    srand(seed)
    n = mat.n
    start = Vector{StartPt}()
    for i in 1:n
        offset = -log(rand())*scale
        push!(start, StartPt(i, offset)) 
    end
    return start
end # makeUniformExpRVStartPoints

function makeDegreeStartPoints{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind}, seed::Int64 = 1)
    srand(seed)
    n = mat.n
    start = Vector{StartPt}()
    for i in 1:n
        wtDeg = 0
        for ind in mat.colptr[i]:(mat.colptr[i+1]-1)
            wt = mat.nzval[ind]
            wtDeg += 1/wt
        end
        offset = -log(rand())*wtDeg
        push!(start, StartPt(i, offset)) 
    end
    return start
end # makeDegreeStartPoints

# Assumes the graph is already connected
# Returns the tree, and the number of edges from each start point
function fractalPrim{Tval,Tind}(mat::SparseMatrixCSC{Tval,Tind}, start::Array{StartPt,1}, seed::Int64 = 1)

    srand(seed)
    n = mat.n
    m = nnz(mat)

    (ai, aj, av) = findnz(mat)

    # flipInd = Laplacians.flipIndex(mat)
    
    # Connected components of graph
    comps = IntDisjointSets(n)

    start_counts = zeros(length(start))
    s = Vector{ValWtPair}()

    treeEdges = zeros(Tind,n-1)
    numEdges = 0
   
    # Add edges
    for i in 1:length(start)
        st = start[i]
        v = st.val
        for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
            offset = st.offset 
            wt = mat.nzval[ind]
            Base.Collections.heappush!(s, ValWtPair(ind, offset-log(rand())/wt, i))
        end
    end

    # More than 1 connected component
    # Will break if graph is not already connected
    while comps.ngroups != 1

      valwt = Base.Collections.heappop!(s) #hog 
      edge = valwt.val

      v = ai[edge]

      if !DataStructures.in_same_set(comps, ai[edge], aj[edge])
          DataStructures.union!(comps, ai[edge], aj[edge])

          numEdges += 1
          treeEdges[numEdges] = edge

          start_counts[valwt.st_pt] += 1
          
          for ind in mat.colptr[v]:(mat.colptr[v+1]-1)
              u = mat.rowval[ind]

              wt = mat.nzval[ind]
              if !DataStructures.in_same_set(comps, v, u)
                  new_pair = 
                  Base.Collections.heappush!(s, ValWtPair(ind, valwt.wt -log(rand())/wt, valwt.st_pt))

              end
          end
      end
    end

    tr = submatrixCSC(mat,treeEdges)
    tr = tr + tr';

    return tr, start_counts

end # fractalPrim

