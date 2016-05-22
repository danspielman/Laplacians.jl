function sampleMany{Tv,Ti}(p::Array{Tv,1}, sampCount::Ti)

    # generate a sorted distribution
    dist = randexp(sampCount + 1)
    dist = dist / sum(dist)
    for i in 1:sampCount
        dist[i + 1] += dist[i]
    end
    pop!(dist)

    n = length(p)
    newp = copy(p) / sum(p)
    for i in 1:(n-1)
        newp[i + 1] += newp[i]
    end

    samples = Array{Ti,1}(sampCount)

    ind = 1
    for i in 1:sampCount
        left = (ind > 1) ? newp[ind - 1] : 0
        right = newp[ind]

        while right < dist[i]
            ind = ind + 1
            left = (ind > 1) ? newp[ind - 1] : 0
            right = newp[ind]
        end

        samples[i] = ind
    end

    return samples
end


#=

    #to test the distribution, we can do something like this

    p = rand(10)
    s = sum(p)

    samp = sampleMany(p, 10000)

    nr = 0
    for i in 1:10000
        if samp[i] == 5
            nr = nr + 1
        end
    end
    
    # then, these two should be very close
    print(nr / 10000)
    print(p[3] / s)

=#