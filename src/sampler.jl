# Sampler:
#
# for keeping a bunch of integers between 1 and n (uses memory n),
# with a positive number for each,
# and sampling with probability proportional to that number
#
# see the notebook for a demo
#
# By Dan on Aug 30, 2015

import Base.push!

mutable struct Sampler
    sums::Array{Float64,1}
    nitems::Int
    n::Int
end

Sampler(n::Int) = Sampler(zeros(Float64,2*n-1),0,n)

# start with this vector
function Sampler(x::Array{Float64,1})
    n = length(x)

    s = Sampler(n)

    s.sums[n:(2*n-1)] = x
    
    for i in (n-1):-1:1
        s.sums[i] = s.sums[2*i] + s.sums[2*i+1]
    end

    s.nitems = n

    return s
end


function push!(s::Sampler, i0::Int, x::Float64)
    i = i0 + s.n - 1

    if s.sums[i] != 0
        error("pushed item " , i , " is already in the sampler.")
    end

    while i > 0
        s.sums[i] = s.sums[i] + x
        i = div(i,2)
    end
    s.nitems = s.nitems + 1
end


function sample(s::Sampler)

    r = rand()*s.sums[1]
    i = 1
    
    # move through sums
    while i < s.n
        x = s.sums[2*i]
        if (r <= x)
            i = 2*i
        else
            i = 2*i+1
            r = r - x
        end
    end

    return i - s.n + 1
end

function remove!(s::Sampler, i::Int)
    i = i + s.n - 1
    x = s.sums[i]
    s.sums[i] = 0
    i = div(i,2)
    while i > 0
        s.sums[i] = s.sums[i] - x
        i = div(i,2)
    end
    s.nitems = s.nitems - 1
end


function pop!(s::Sampler)
    if s.nitems < 1
        error("No items to pop")
    end
    
    i = sample(s)
    remove!(s,i)
    return i
end

