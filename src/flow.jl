" implementation of Dinic's algorithm. computes the maximum flow and min-cut in G between s and t. we consider the adjacency matrix to be the capacity matrix "
function maxflow{Tv, Ti}(G::SparseMatrixCSC{Tv, Ti}, s::Int64, t::Int64; justflow = true)

  n = max(G.n, G.m)

  # initialize the capacity matrix
  C = copy(G)
  BI = backIndices(G)

  Q = zeros(Ti, n)
  inQ = zeros(Bool, n)
  prev = [(0,0) for i in 1:n]
  totalflow = 0

  improve = true
  while improve
    improve = false

    # initialize the queue at the current step
    left = 0
    right = 1
    Q[right] = s
    inQ[s] = true

    # do bfs
    while left < right
      left = left + 1
      u = Q[left]
      if (u == t)
        continue
      end

      for i in 1:deg(G, u)
        v = nbri(G, u, i)

        if inQ[v] == false && weighti(C, u, i) > 0
          right = right + 1
          Q[right] = v
          inQ[v] = true
          prev[v] = (u, i)
        end
      end
    end

    # do block flow
    if inQ[t] == true
      improve = true

      for i in 1:deg(G, t)
        u = nbri(G, t, i)

        currentflow = C[u, t]
        while u != s
          prevU = prev[u][1]
          indexPrevU = prev[u][2]
          currentflow = min(currentflow, weighti(BI, prevU, indexPrevU))
          u = prevU
        end

        if currentflow != 0
          # let's consider that prev[t] = (u, i)
          u = nbri(G, t, i)
          prev[t] = (u, weighti(BI, t, i))
          u = t

          #=
            prevU is the vertex immmediately before u on the path from s to t
            indexU is such that nbri(G, u, indexU) = prevU
            indexPrevU is such that nbri(G, prevU, indexPrevU) = u
            capPrevUU and capUPrevU have similar meanings
          =#
          while u != s
            prevU = prev[u][1]
            indexPrevU = prev[u][2]
            indexU = weighti(BI, prevU, indexPrevU)

            capPrevUU = weighti(C, prevU, indexPrevU)
            capUPrevU = weighti(C, u, indexU)

            setValue(C, prevU, indexPrevU, capPrevUU - currentflow)
            setValue(C, u, indexU, capUPrevU + currentflow)

            u = prevU
          end
        end

        totalflow = totalflow + currentflow
      end

    end

    # reset the queue
    for i in 1:right
      u = Q[i]

      Q[i] = 0
      inQ[u] = false
      prev[u] = (0,0)
    end
    left = right = 0

  end

  if justflow == true
    return totalflow
  else
    # compute the min-cut
    left = 0
    right = 1
    Q[right] = s
    inQ[s] = true

    while left < right
      left = left + 1

      u = Q[left]
      for i in 1:deg(G, u)
        v = nbri(G, u, i)

        if inQ[v] == 0 && weighti(C, u, i) > 0
          right = right + 1
          Q[right] = v
          inQ[v] = 1
        end
      end
    end

    cut = Tuple{Int64,Int64}[]
    for (i in 1:right)
      u = Q[i]
      for j in 1:deg(G, u)
        v = nbri(G, u, j)

        if weighti(C, u, j) == 0
          push!(cut, (u, v))
        end
      end
    end

    return totalflow, cut
  end

end
