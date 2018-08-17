"""
    a = rand_regular_bipartite(n,k)

Random k-regular bipartite graph between two sets of n vertices.
No repeat edges, so can take a long time to build of k is close to n.

Returns a (possibly) asymmetric matrix that contains the upper-right block.
"""
function rand_regular_bipartite(n,k)
    @assert k <= div(n,2)
    
    p = randperm(n)
    
    # start by building the upper-right block    
    b = sparse(1:n, p, 1, n, n) 
    
    for i in 2:k
        done = false
        while !done
            p = randperm(n)
            done = true
            for i in 1:n
                if b[i,p[i]] > 0
                    done = false
                end
            end
            if done 
                b = b + sparse(1:n, p, 1, n, n)
            end
        end
        
    end
    
    return b    
end

function test_regular_bipartite(B)
    d = maximum(sum(B,dims=1))
    d == minimum(sum(B,dims=1)) &&
    d == minimum(sum(B,dims=2)) &&
    d == maximum(sum(B,dims=2)) &&
    length(unique(B)) <= 2
end

function regular_bipartite_matching(B)
    @assert test_regular_bipartite(B)
    
    n = size(B,1)
    match = zeros(Int,n) 
    
    # match is indexed by right vertices, and gives name of match on left
    
    matched = zeros(Bool,n) # indicates which vertices on left are matched
    
    d = maximum(sum(B,dims=1))    
    if d == 1
        return invperm(B.rowval)
    end
    
    for i in 1:n
        ops = match_walk(B, match, matched, n, d)
        
        # eliminate loops in ops
        us = [op[1] for op in ops]

        firsts = indexin(us,us)
        i = length(ops)
        while i > 0
            (u,v) = ops[i]
            match[v] = u
            matched[u] = true
            
            if firsts[i] < i
                i = firsts[i] - 1
            else
                i = i - 1
            end
        end
        
    end
    
    return match
end

function match_walk(B, match, matched, n, d)
    
    ops = Tuple{Int64,Int64}[]
   
    u = rand(findall(.!matched))
    v = 0

    
    done = false
    while !done
        
        foundv = false
        while !foundv
            i = rand(1:d)
            v = B.rowval[B.colptr[u]+i-1]
            foundv = (u != match[v])
        end
        
        push!(ops,(u,v))
        
        if match[v] == 0
            done = true
        else
            u = match[v]
        end
        
    end
    
    return ops
end

"""
    S = latin_square(n)

Computes a random (but not uniformly random) n-by-n latin square.
"""
function latin_square(n)
    
    S = zeros(Int,n,n)
    
    B = ones(Int,n,n)
    
    for j in 1:n
   
        match = regular_bipartite_matching(sparse(B))
        for i in 1:n
            S[i,match[i]] = j
            B[i,match[i]] = 0
        end        
    end
    
    return S
end

function test_latin_square(S)
    n = size(S,1)
    
    b = true
    b = b && size(S,2) == n
    
    for i in 1:n
        b = b && length(unique(S[:,i])) == n
        b = b && length(unique(S[i,:])) == n
    end
    
    b = b && length(unique(S)) == n
    return b
end
    

