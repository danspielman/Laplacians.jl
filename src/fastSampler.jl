# Alias sampling  
immutable FastSampler{Tv,Ti}
    F::Array{Tv,1}    # array of p[s.A[i]] / (1 / n)
    A::Array{Ti,1}    # array of indices of elements < (1 / n)
    V::Array{Ti,1}    # arary of indices of elements > (1 / n)
    n::Ti
    rng::AbstractRNG  # the pseudo random number generator
end

# sample a random number
function sample{Tv,Ti}(s::FastSampler{Tv,Ti})
    #i = rand(1:s.n)
    i = ceil(Ti,rand(s.rng)*s.n)
    f = rand(s.rng)
    if f < s.F[i]
        return s.A[i]
    else
        return s.V[i]
    end
end

# sample a random number
function sampleMany{Tv,Ti}(s::FastSampler{Tv,Ti},sampCount::Ti)
    #i = rand(1:s.n)
    samples = Array{Ti,1}(sampCount)
    for j = 1:sampCount
        i = ceil(Ti,rand(s.rng)*s.n)
        f = rand(s.rng)
        if f < s.F[i]
            samples[j] = s.A[i]
        else
            samples[j] = s.V[i]
        end
    end
    return samples
end

"""
    s = FastSampler(p)
"""
function FastSampler{Tv}(p::Array{Tv,1}; residual::Bool = false, rng::AbstractRNG=Base.Random.GLOBAL_RNG)

    @assert(minimum(p) > 0, "The probability vector has a negative entry")

    n = length(p)

    # generate the Sampler
    F = Array{Tv}(n)
    A = Array{Int64}(n)
    V = Array{Int64}(n)
    
    # generate a normalized p, so that it sums up to n
    newp = n * copy(p) / sum(p)

    # fill up the arays A and V
    pos = 0                      # number of filled or half filled buckets
    posFilled = 0                # number of filled buckets

    for i in 1:n
        if newp[i] < 1 #&& !epsequal(newp[i], 1)
            pos = pos + 1
            A[pos] = i
            F[pos] = newp[i]
        end
    end

    err = 0
    for i in 1:n
        if newp[i] >= 1 
            val = newp[i]
            
            while val >= 1 
                if posFilled < pos
                    # fill up a bucket already containing a small element, and subtract that value from the current element
                    posFilled = posFilled + 1
                    V[posFilled] = i
                    val = val - (1 - F[posFilled])
                else
                    # create a bucket with only one large element. note posFilled = pos
                    pos = pos + 1
                    posFilled = posFilled + 1
                    A[pos] = V[pos] = i
                    F[pos] = 1
                    val = val - 1
                end
            end

            if val > 0 && pos == n
            	err = err + val
            end

            # treat the case in which val becomes smaller than 1. Have to do the pos < n check because of precision errors
            if val > 0 && pos < n
                pos = pos + 1
                A[pos] = i
                F[pos] = val
            end
        end
    end

    while posFilled < n
        # if we are in this case, posFilled should be equal to pos - 1
        posFilled = posFilled + 1
        V[posFilled] = A[posFilled]
        err += 1-F[posFilled]
    end

    if residual
        return FastSampler(F, A, V, n, rng), err
    else
        return FastSampler(F, A, V, n, rng)
    end
end


"""
    s = blockSample(p; k = length(p))

Compute numbers sampled with probability proportional to the vector p.
They are returned in order.
If need be, they can be permuted with randperm.
"""
function blockSample(p; k = length(p))
    samp = zeros(Int,k)

    r = sort(rand(k))

    cs = cumsum(p)
    cs = cs / cs[end]
    
    j = 1
    for i in 1:k
        while r[i] > cs[j]
            j += 1
        end

        samp[i] = j
    end

    return samp
end
