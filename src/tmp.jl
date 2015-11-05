function randWeight(a)

    n = a.n
    (ai,aj) = findnz(a)
    m = length(ai)
    
    # potentials or edge-based

    if (rand(1)[1] < .3)
        w = rand(m)

    else
        v = randn(a.n)

        # mult by matrix ?
        if (rand(1)[1] < .5)

            invdeg = spdiagm(1./(a*ones(size(a)[1])))
            if (rand(1)[1] < .5)
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

    if (rand(1)[1] < .5)
        w = 1./w
    end

    w = w / mean(w)

    ar = sparse(ai,aj,w,n,n)
    ar = ar + ar';
    return ar
end