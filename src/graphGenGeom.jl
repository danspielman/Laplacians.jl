using Random
using Polynomials

# mutable struct IJV3d{Tv,Ti}
#     n::Ti
#     nnz::Ti
#     i::Array{Tuple{Ti,Ti,Ti},1}
#     j::Array{Tuple{Ti,Ti,Ti},1}
#     v::Array{Tv,1}
# end

mutable struct IJV3d{Tv,Ti}
    nnz::Ti
    i::Array{Tuple{Ti,Ti,Ti},1}
    j::Array{Tuple{Ti,Ti,Ti},1}
    v::Array{Tv,1}
    xlen::Ti
    ylen::Ti
end

function indexToLinear(x, y, xlen)
    return x + (y - 1) * xlen
end

function indexToLinear(x, y, z, xlen, ylen)
    return x + (y - 1) * xlen + (z - 1) * xlen * ylen
end

function indexToLinear(p::Tuple{Ti,Ti,Ti}, xlen::Integer, ylen::Integer) where {Ti}
    return p[1] + (p[2] - 1) * xlen + (p[3] - 1) * xlen * ylen
end

function ggrid2_ijv(s1::Integer, s2::Integer)
    n = s1 * s2
    nnz = 2 * n - s1 - s2
    
    i = Array{Int64,1}(undef, nnz)
    j = Array{Int64,1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1

    for y = 1:s2
        for x = 1:s1
            # horizontal
            if x < s1
                i[e] = indexToLinear(x, y, s1)
                j[e] = indexToLinear(x + 1, y, s1)
                v[e] = 1
                e += 1  
            end 
            # vertical
            if y < s2
                i[e] = indexToLinear(x, y, s1)
                j[e] = indexToLinear(x, y + 1, s1)
                v[e] = 1
                e += 1
            end 
        end
    end
    return IJV(n, nnz, i, j, v)
end

function ggrid2(s1::Integer, s2::Integer)
    A = sparse(ggrid2_ijv(s1, s2))
    return A + A'
end

function ggrid2_checkered_ijv(s1::Integer, s2::Integer, b1::Integer, b2::Integer, w::Real)
    l1 = div(s1, b1)
    l2 = div(s2, b2)

    n = s1 * s2
    nnz = 2 * n - s1 - s2
    
    i = Array{Int64,1}(undef, nnz)
    j = Array{Int64,1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1

    for x = 1:s1
        for y = 1:s2
            checkered = (mod(div(x - 1, l1) + div(y - 1, l2), 2) == 0)
            # horizontal
            if x < s1
                i[e] = indexToLinear(x, y, s1)
                j[e] = indexToLinear(x + 1, y, s1)
                v[e] = checkered ? w : 1
                e += 1
            end 
            # vertical
            if y < s2
                i[e] = indexToLinear(x, y, s1)
                j[e] = indexToLinear(x, y + 1, s1)
                v[e] = checkered ? w : 1
                e += 1
            end 
        end
    end
    return IJV(n, nnz, i, j, v)
end

function ggrid2_checkered(s1::Integer, s2::Integer, b1::Integer, b2::Integer, w::Real)
    A = sparse(ggrid2_checkered_ijv(s1, s2, b1, b2, w))
    return A + A'
end

function ggrid2_checkered(s::Integer, b::Integer, w::Real)
    return ggrid2_checkered(s, s, b, b, w)
end

function plot_graph_weighted(gr, x, y;dots = true,setaxis = true,number = false)

    if isa(color, Vector) && length(color) == 3
        col = RGB(color...)
    else
        col = color
    end

    p = plot(;legend = false, axis = false, xticks = false, yticks = false)

    (ai, aj, av) = findnz(triu(gr))
    for i in 1:length(ai)
        s = [ai[i]; aj[i]]
        plot!(p, x[s], y[s], line_z = av[i], linecolor = :blues)
    end

    if dots
        scatter!(p, x, y, markercolor = col, markerstrokecolor = false)
    end

    if number
        annotate!(p, x, y, collect(1:length(x)))
    end

    minx = minimum(x)
    maxx = maximum(x)
    miny = minimum(y)
    maxy = maximum(y)
    delx = maxx - minx
    dely = maxy - miny

    plot!(p; ylims = (miny - dely / 20, maxy + dely / 20))
    plot!(p; xlims = (minx - delx / 20, maxx + delx / 20))

    display(p)
    return p
end

function ggrid2coords(s1::Integer, s2::Integer)
    x = kron(ones(s2), collect(1:s1))
    y = kron(collect(1:s2), ones(s1))
    return x, y
end # grid2coords


function ggrid3_ijv(s1::Integer, s2::Integer, s3::Integer)
    n = s1 * s2 * s3
    nnz = 3 * n - s1 * s2 - s1 * s3 - s2 * s3
    
    i = Array{Int64,1}(undef, nnz)
    j = Array{Int64,1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1

    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                # eastward edge
                if x < s1
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x + 1, y, z, s1, s2)
                    v[e] = 1
                    e += 1  
                end 
                # northward edge
                if y < s2
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y + 1, z, s1, s2)
                    v[e] = 1
                    e += 1
                end 
                # upward edge
                if z < s3
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y, z + 1, s1, s2)
                    v[e] = 1
                    e += 1
                end 
            end
        end
    end
    return IJV(n, nnz, i, j, v)
end

function ggrid3_ijv3d(s1::Integer, s2::Integer, s3::Integer)
    n = s1 * s2 * s3
    nnz = 3 * n - s1 * s2 - s1 * s3 - s2 * s3
    
    i = Array{Tuple{Int64,Int64,Int64},1}(undef, nnz)
    j = Array{Tuple{Int64,Int64,Int64},1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1

    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                # eastward edge
                if x < s1
                    i[e] = (x, y, z) #indexToLinear(x,y,z,s1,s2)
                    j[e] = (x + 1, y, z) #indexToLinear(x+1,y,z,s1,s2)
                    v[e] = 1
                    e += 1  
                end 
                # northward edge
                if y < s2
                    i[e] = (x, y, z) #indexToLinear(x,y,z,s1,s2)
                    j[e] = (x, y + 1, z) #indexToLinear(x,y+1,z,s1,s2)
                    v[e] = 1
                    e += 1
                end 
                # upward edge
                if z < s3
                    i[e] = (x, y, z) #indexToLinear(x,y,z,s1,s2)
                    j[e] = (x, y, z + 1) #indexToLinear(x,y,z+1,s1,s2)
                    v[e] = 1
                    e += 1
                end 
            end
        end
    end
    return IJV3d(nnz, i, j, v, s1, s2)
end

function halfopenggrid3_ijv3d(s1::Integer, s2::Integer, s3::Integer)
    n = s1 * s2 * s3
    nnz = 3 * n - s1 * s2 - s1 * s3 - s2 * s3
    
    i = Array{Tuple{Int64,Int64,Int64},1}(undef, nnz)
    j = Array{Tuple{Int64,Int64,Int64},1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1
    # TODO SOMETHING WRONG? -- maybe with number of edges allocated?
    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                if x < s1 && y < s2 && z < s3
                    # eastward edge
                    i[e] = (x, y, z) 
                    j[e] = (x + 1, y, z) 
                    v[e] = 1
                    e += 1  

                     # northward edge
                     i[e] = (x, y, z) 
                     j[e] = (x, y + 1, z) 
                     v[e] = 1
                     e += 1

                     # upward edge
                     i[e] = (x, y, z) 
                     j[e] = (x, y, z + 1) 
                     v[e] = 1
                     e += 1
                end 
            end
        end
    end
    return IJV3d(nnz, i, j, v, s1, s2)
end

function shift(a::IJV3d{Tv,Ti}, s1::Ti, s2::Ti, s3::Ti) where {Tv,Ti}
    i_shifted = map(x->(x[1] + s1, x[2] + s2, x[3] + s3), a.i)
    j_shifted = map(x->(x[1] + s1, x[2] + s2, x[3] + s3), a.j)
    return IJV3d(a.nnz, i_shifted, j_shifted, a.v, a.xlen + s1, a.ylen + s2)
end

function extend(a::IJV3d{Tv,Ti}, s1::Ti, s2::Ti, s3::Ti) where {Tv,Ti}
    IJV3d(a.nnz, a.i, a.j, a.v, a.xlen + s1, a.ylen + s2)
end

function add_graphs_3d(a::IJV3d{Tv,Ti}, b::IJV3d{Tv,Ti}) where {Tv,Ti}
    @assert(a.xlen == b.xlen)
    @assert(a.ylen == b.ylen)
    return IJV3d(a.nnz + b.nnz, [a.i; b.i], [a.j; b.j], [a.v; b.v], a.xlen, a.ylen) 
end

# this is just to test splitting
function ggrid3_testsplit(s1::Integer, s2::Integer, s3::Integer)
    s1a = div(s1 + 1, 2)
    s1b = s1+1 - s1a
    a = halfopenggrid3_ijv3d(s1a, s2, s3) 
    b = halfopenggrid3_ijv3d(s1b, s2, s3)
    amod = extend(a,s1b-1,0,0);
    bmod = shift(b,s1a-1,0,0);
    return add_graphs_3d(amod, bmod)
    #note: the grids need to share the  middle vertex
end

function ijv3d_to_sparse(a::IJV3d{Tv,Ti})where {Tv, Ti}
    i_shifted = map(p -> indexToLinear(p, a.xlen, a.ylen),a.i)
    j_shifted = map(p -> indexToLinear(p, a.xlen, a.ylen),a.j)
    return sparse(i_shifted, j_shifted, a.v)
end

function ggrid3_checkered_ijv(s1::Integer, s2::Integer, s3::Integer, b1::Integer, b2::Integer, b3::Integer, w::Real)
    l1 = div(s1, b1)
    l2 = div(s2, b2)
    l3 = div(s3, b3)

    n = s1 * s2 * s3
    nnz = 3 * n - s1 * s2 - s1 * s3 - s2 * s3
    
    i = Array{Int64,1}(undef, nnz)
    j = Array{Int64,1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1

    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                checkered = (mod(div(x - 1, l1) + div(y - 1, l2) + div(z - 1, l3), 2) == 0)

                # eastward edge
                if x < s1
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x + 1, y, z, s1, s2)
                    v[e] = checkered ? w : 1
                    e += 1  
                end 
                # northward edge
                if y < s2
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y + 1, z, s1, s2)
                    v[e] = checkered ? w : 1
                    e += 1
                end 
                # upward edge
                if z < s3
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y, z + 1, s1, s2)
                    v[e] = checkered ? w : 1
                    e += 1
                end 
            end
        end
    end
    return IJV(n, nnz, i, j, v)
end

function ggrid3_aniso_ijv(s1::Integer, s2::Integer, s3::Integer, w1::Real, w2::Real=1, w3::Real=1)
    n = s1 * s2 * s3
    nnz = 3 * n - s1 * s2 - s1 * s3 - s2 * s3
    
    i = Array{Int64,1}(undef, nnz)
    j = Array{Int64,1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1

    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                if x < s1
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x + 1, y, z, s1, s2)
                    v[e] = w1
                    e += 1  
                end 
                # northward edge
                if y < s2
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y + 1, z, s1, s2)
                    v[e] = w2
                    e += 1
                end 
                # upward edge
                if z < s3
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y, z + 1, s1, s2)
                    v[e] = w3
                    e += 1
                end 
            end
        end
    end
    return IJV(n, nnz, i, j, v)
end

function ggrid3(s1::Integer, s2::Integer, s3::Integer)
    A = sparse(ggrid3_ijv(s1, s2, s3))
    return A + A'
end


function ggrid3_aniso(s1::Integer, s2::Integer, s3::Integer, w1::Real, w2::Real=1, w3::Real=1)
    A = sparse(ggrid3_aniso_ijv(s1, s2, s3, w1, w2, w3))
    return A + A'
end

function ggrid3_checkered(s1::Integer, s2::Integer, s3::Integer, b1::Integer, b2::Integer, b3::Integer, w::Real)
    A = sparse(ggrid3_checkered_ijv(s1, s2, s3, b1, b2, b3, w))
    return A + A'
end

function ggrid3coords(s1::Integer, s2::Integer, s3::Integer)
    x = kron(ones(s3), kron(ones(s2), collect(1:s1)))
    y = kron(ones(s3), kron(collect(1:s2), ones(s1)))
    z = kron(collect(1:s3), kron(ones(s2), ones(s1)))
    return x, y, z
end # grid2coords

function plot_graph_weighted(gr, x, y, z)#;dots=true,setaxis=true)#,number=false)

    if isa(color, Vector) && length(color) == 3
        col = RGB(color...)
    else
        col = color
    end

    p = plot3d(;legend = false, axis = false, xticks = false, yticks = false, zticks = false)

    (ai, aj, av) = findnz(triu(gr))
    for i in 1:length(ai)
        s = [ai[i]; aj[i]]
        plot3d!(p, x[s], y[s], z[s], line_z = av[i], linecolor = :blues)
    end

    # if dots
    #     scatter!(p, x, y, z, markercolor=col, markerstrokecolor=false)
    # end

    # if number
    #     annotate3d!(p, x, y, z, collect(1:length(x)))
    # end
    # 3d annotations currently not supported in Plots, but maybe in Plotly?
    # https://stackoverflow.com/questions/41333965/how-to-make-3d-annotations-in-julia-plots
    # https://stackoverflow.com/questions/37596865/add-annotation-to-3d-scatterplot-in-plotly

    minx = minimum(x)
    maxx = maximum(x)
    miny = minimum(y)
    maxy = maximum(y)
    delx = maxx - minx
    dely = maxy - miny

    # RK: what did this do before?
    #plot!(p; ylims = (miny - dely/20, maxy + dely/20))
    #plot!(p; xlims = (minx - delx/20, maxx + delx/20))

    display(p)
    return p
end

function getBoundary2(s1::Integer, s2::Integer)
    n = s1 * s2
    bot = 1:s1
    top = n - s1 + 1:n
    left = s1 + 1:s1:n - s1
    right = 2 * s1:s1:n - 1
    bndry = [bot; top; left; right]
    return bndry
end

function getInterior2(s1::Integer, s2::Integer)
    n = s1 * s2
    return setdiff(1:n, getBoundary2(s1, s2))
end

function getBoundary3(s1::Integer, s2::Integer, s3::Integer)
    n = s1 * s2 * s3
    zsidesize = s1 * s2
    xsidesize = s2 * s3
    ysidesize = s1 * s3

    # alternatively, we could get z boundary as
    # zside1b = 1:zsidesize # bot 
    # zside2b = n-zsidesize+1:n # top
    
    zside1 = Array{Int64}(undef, zsidesize)
    for (i, (x, y)) in Base.Iterators.enumerate(Base.product(1:s1, 1:s2))
        zside1[i] = indexToLinear(x, y, 1, s1, s2)
    end
    
    zside2 = Array{Int64}(undef, zsidesize)
    for (i, (x, y)) in Base.Iterators.enumerate(Base.product(1:s1, 1:s2))
        zside2[i] = indexToLinear(x, y, s3, s1, s2)
    end
    
    xside1 = Array{Int64}(undef, xsidesize)
    for (i, (y, z)) in Base.Iterators.enumerate(Base.product(1:s2, 1:s3))
        xside1[i] = indexToLinear(1, y, z, s1, s2)
    end
    
    xside2 = Array{Int64}(undef, xsidesize)
    for (i, (y, z)) in Base.Iterators.enumerate(Base.product(1:s2, 1:s3))
        xside2[i] = indexToLinear(s1, y, z, s1, s2)
    end
    
    yside1 = Array{Int64}(undef, ysidesize)
    for (i, (x, z)) in Base.Iterators.enumerate(Base.product(1:s1, 1:s3))
        yside1[i] = indexToLinear(x, 1, z, s1, s2)
    end
    
    yside2 = Array{Int64}(undef, ysidesize)
    for (i, (x, z)) in Base.Iterators.enumerate(Base.product(1:s1, 1:s3))
        yside2[i] = indexToLinear(x, s2, z, s1, s2)
    end
    
    return [zside1; zside2; xside1; xside2; yside1; yside2]
end

function getInterior3(s1::Integer, s2::Integer, s3::Integer)
    n = s1 * s2 * s3
    return setdiff(1:n, getBoundary3(s1, s2, s3))
end

# function zeroBoundaryLap2(s1,)



function ggrid3_checkeredrandom_ijv(s1::Integer, s2::Integer, s3::Integer, b1::Integer, b2::Integer, b3::Integer, multrange::Integer, seed::Integer = 1234)
    l1 = div(s1, b1)
    l2 = div(s2, b2)
    l3 = div(s3, b3)
    b1odd = mod(b1, 2) == 0 ? b1 + 1 : b1
    b2odd = mod(b2, 2) == 0 ? b2 + 1 : b2
    b3odd = mod(b3, 2) == 0 ? b3 + 1 : b3

    n = s1 * s2 * s3
    nnz = 3 * n - s1 * s2 - s1 * s3 - s2 * s3
    
    i = Array{Int64,1}(undef, nnz)
    j = Array{Int64,1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1

    # rng = MersenneTwister(seed);

    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                boxInd1 = div(x - 1, l1) + 1
                boxInd2 = div(y - 1, l2) + 1
                boxInd3 = div(z - 1, l3) + 1
                boxIndLinear = indexToLinear(boxInd1, boxInd2, boxInd3, b1odd, b2odd)
                boxIndLinearCheck = mod(boxIndLinear - 1, 2) == 0
                checkered = (mod(div(x - 1, l1) + div(y - 1, l2) + div(z - 1, l3), 2) == 0)
                @assert(boxIndLinearCheck == checkered)
                wexponent = mod(hash(seed + boxIndLinear), multrange)
                w = 10.0^(-wexponent)

                # eastward edge
                if x < s1
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x + 1, y, z, s1, s2)
                    v[e] = checkered ? w : 1
                    e += 1  
                end 
                # northward edge
                if y < s2
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y + 1, z, s1, s2)
                    v[e] = checkered ? w : 1
                    e += 1
                end 
                # upward edge
                if z < s3
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y, z + 1, s1, s2)
                    v[e] = checkered ? w : 1
                    e += 1
                end 
            end
        end
    end
    return IJV(n, nnz, i, j, v)
end

function ggrid3_checkeredrandom(s1::Integer, s2::Integer, s3::Integer, b1::Integer, b2::Integer, b3::Integer, multrange::Integer, seed::Integer = 1234)
    A = sparse(ggrid3_checkeredrandom_ijv(s1, s2, s3, b1, b2, b3, multrange, seed))
    return A + A'
end

function ggrid3_perc_ijv(s1::Integer, s2::Integer, s3::Integer, w::Real, seed::Integer = 1234, p::Real = 0.2488)
    n = s1 * s2 * s3
    nnz = 3 * n - s1 * s2 - s1 * s3 - s2 * s3
    
    i = Array{Int64,1}(undef, nnz)
    j = Array{Int64,1}(undef, nnz)
    v = Array{Float64,1}(undef, nnz)

    e = 1

    rng = MersenneTwister(seed);

    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                if x < s1
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x + 1, y, z, s1, s2)
                    v[e] = rand(rng) < p ? w : 1
                    e += 1  
                end 
                # northward edge
                if y < s2
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y + 1, z, s1, s2)
                    v[e] = rand(rng) < p ? w : 1
                    e += 1
                end 
                # upward edge
                if z < s3
                    i[e] = indexToLinear(x, y, z, s1, s2)
                    j[e] = indexToLinear(x, y, z + 1, s1, s2)
                    v[e] = rand(rng) < p ? w : 1
                    e += 1
                end 
            end
        end
    end
    return IJV(n, nnz, i, j, v)
end


function ggrid3_perc(s1::Integer, s2::Integer, s3::Integer, w::Real, seed::Integer = 1234, p::Real = 0.2488)
    A = sparse(ggrid3_perc_ijv(s1, s2, s3, w, seed, p))
    return A + A'
end

function find_nx(nnz::Integer, xi::Real)
    p = Polynomial([-nnz, 0, - 2 * (xi+2)/xi, 7/xi])
    r = roots(p)
    for root in r
        if imag(root) == 0.0
            return real(root)
        end
    end
end


"""
Wrapper for uniform grid.
It returns the sddm of uniform grid with the number of entries equals nnz.
"""
function uniform_grid_sddm(nnz::Integer) 
    frac_nx = find_nx(nnz, 1)
    nx = round(Integer, frac_nx)
    nz = round(Integer, frac_nx)
    A = ggrid3_aniso(nx+2,nx+2,nz+2,1)
    L = lap(A)
    int = getInterior3(nx+2,nx+2,nz+2)
    M = L[int,int]
    return M 
end


"""
Wrapper for checkered grid.
It returns the sddm of checkered grid with the number of entries equals nnz.
"""
function checkered_grid_sddm(nnz::Integer, b1::Integer, b2::Integer, b3::Integer, w::Real)
    frac_nx = find_nx(nnz, 1)
    s = round(Integer, frac_nx) + 2
    A = ggrid3_checkered(s, s, s, b1, b2, b3, w)
    L = lap(A)
    int = getInterior3(s, s, s)
    M = L[int, int]
    return M
end

"""
Wrapper for anisotropic grid.
It returns the sddm of anisotropic grid with the number of entries equals nnz.
"""

function aniso_grid_sddm(nnz::Integer, xi::Real)
    frac_nx = find_nx(nnz, xi)
    nx = round(Int, frac_nx)
    nz = round(Int, frac_nx/xi)
    A = ggrid3_aniso(nx+2,nx+2,nz+2,1)
    L = lap(A)
    int = getInterior3(nx+2,nx+2,nz+2)
    M = L[int,int]
    return M
end

"""
Wrapper for weighted grid
It returns the sddm of weighted grid with the number of entries eqauls nnz
"""

function wgrid_sddm(nnz::Integer, w::Real)
    frac_nx = find_nx(nnz, 1)
    nx = round(Int, frac_nx)
    A = 1/w * ggrid3_aniso(nx+2,nx+2,nx+2,w)
    L = lap(A)
    int = getInterior3(nx+2,nx+2,nx+2)
    M = L[int,int]
    return M
end

