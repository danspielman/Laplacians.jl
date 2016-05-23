#TODO: would be nice to add a link to the paper rasmus was referencing    
immutable NewSampler{Tv,Ti}
    F::Array{Tv,1}    # array of p[s.A[i]] / (1 / n)
    A::Array{Ti,1}    # array of indices of elements < (1 / n)
    V::Array{Ti,1}    # arary of indices of elements > (1 / n)
    n::Ti
end

# immutable Sampler
#     F::Array{Float64,1}    # array of p[s.A[i]] / (1 / n)
#     A::Array{Int64,1}    # array of indices of elements < (1 / n)
#     V::Array{Int64,1}    # arary of indices of elements > (1 / n)
#     n
# end

# sample a random number
function newSample{Tv,Ti}(s::NewSampler{Tv,Ti})
    #i = rand(1:s.n)
    i = ceil(Ti,rand()*s.n)
    f = rand()
    if f < s.F[i]
        return s.A[i]
    else
        return s.V[i]
    end

end

# sample a random number
# 
function newSampleMany{Tv,Ti}(s::NewSampler{Tv,Ti},sampCount::Ti)
    #i = rand(1:s.n)
    samples = Array{Ti,1}(sampCount)
    for j = 1:sampCount
        i = ceil(Ti,rand()*s.n)
        f = rand()
        if f < s.F[i]
            samples[j] = s.A[i]
        else
            samples[j] = s.V[i]
        end
    end
    return samples
end

# initialize the sampler. To get the residual error after building the sampler, set residual to true
function newSampler{Tv}(p::Array{Tv,1}; residual::Bool = false)

    n = length(p)

    # generate the Sampler
    F = Array{Tv}(n)
    A = Array{Int64}(n)
    V = Array{Int64}(n)
    
    # generate a normalized p, so that it sums up to n
    newp = n * copy(p) / sum(p)

    # fill up the arays A and V
    pos = 0                        # number of filled or half filled buckets
    posFilled = 0                # number of filled buckets

    for i in 1:n
        if newp[i] < 1 #&& !epsequal(newp[i], 1)
            pos = pos + 1
            A[pos] = i
            F[pos] = newp[i]
        end
    end

    #residualError = 0

    # printErr = 0
    
    for i in 1:n
        if newp[i] >= 1 #|| epsequal(newp[i], 1)
            val = newp[i]

            # if pos == n && i < n
            #     printErr = 1
            # end
            
            while val >= 1 #|| epsequal(val, 1)
                if posFilled < pos
                    # fill up a bucket already containing a small element, and subtract that value from the current element
                    posFilled = posFilled + 1
                    V[posFilled] = i
                    val = val - (1 - F[posFilled])
                else
                    # if pos == n
                    #     writedlm("boundsErr.txt",p)
                    # end
                    
                    # create a bucket with only one large element. note posFilled = pos
                    pos = pos + 1
                    posFilled = posFilled + 1
                    A[pos] = V[pos] = i
                    F[pos] = 1
                    val = val - 1
                end
            end

            # treat the case in which val becomes smaller than 1. Have to do the pos < n check because of precision errors
            if val > 0 && pos < n #&& !epsequal(val, 0) && pos < n
                pos = pos + 1
                A[pos] = i
                F[pos] = val
            end

            # residualError += val
        end
    end

    err = 0
    while posFilled < n
        # if we are in this case, posFilled should be equal to pos - 1
        posFilled = posFilled + 1
        V[posFilled] = A[posFilled]
        err += 1-F[posFilled]
    end

    # if printErr != 0
    #     println("residual error = ", err)
    # end

    #println("err = ",err)
    
    #@assert(posFilled == pos, "pos and posFilled differ")

    if residual
        return NewSampler(F, A, V, n), err
    else
        return NewSampler(F, A, V, n)
    end
end

function epsequal(a, b; eps=1e-17)
    return abs(a - b) <= eps
end
