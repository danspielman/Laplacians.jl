#=
Written by Rasmus Kyng.

Given an adjacency matrix mat
where A[i,j] != 0 indicates a directed edge from i to j,
compute a topological sorting of the vertices, s.t. 

Please check that the edge direction agrees with what you want before usage!
=#

function toposort{Tv,Ti}(mat::SparseMatrixCSC{Tv,Ti})
    # each COL i gives the outgoing edges of vertex i
    n = mat.n
    colptr::Array{Ti,1} = mat.colptr
    rowval::Array{Ti,1} = mat.rowval # the row index corresponding to a given colptr entry
    popSorted::Array{Ti,1} = zeros(Ti,n);
    seen::Array{Bool,1} = falses(n);
    popCount::Ti = 0
    s = DataStructures.Stack(Tuple{Ti,Ti})
    for i = 1:n
        if seen[i]
            continue
        end
        seen[i] = true
        DataStructures.push!(s,(i,colptr[i]))
        while !isempty(s)
            vertNbrTuple = DataStructures.pop!(s)
            v = vertNbrTuple[1]
            edgeNum = vertNbrTuple[2]
            if edgeNum < colptr[v+1]
                DataStructures.push!(s,(v,edgeNum+1))
                nbr = rowval[edgeNum]
                if !seen[nbr]
                    seen[nbr] = true
                    DataStructures.push!(s,(nbr,colptr[nbr]))
                end
            else
                popCount += 1
                popSorted[popCount] = v
            end
        end
    end
    return popSorted
end


"""The signed edge-vertex adjacency matrix"""
function dirEdgeVertexMat(A::SparseMatrixCSC)
    (ai,aj) = findnz(A)
    m = length(ai)
    n = size(A)[1]
    return sparse(collect(1:m),ai,1.0,m,n) - sparse(collect(1:m),aj,1.0,m,n)
end

