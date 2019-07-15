function indexToLinear(x,y,xlen)
    return x+(y-1)*xlen
end

function indexToLinear(x,y,z,xlen,ylen)
    return x+(y-1)*xlen+(z-1)*xlen*ylen
end

function ggrid2_ijv(s1::Integer,s2::Integer)
    n = s1*s2
    nnz = 2*n - s1 - s2
    
    i = Array{Int64,1}(undef,nnz)
    j = Array{Int64,1}(undef,nnz)
    v = Array{Float64,1}(undef,nnz)

    e = 1

    for y = 1:s2
        for x = 1:s1
            # horizontal
            if x < s1
                i[e] = indexToLinear(x,y,s1)
                j[e] = indexToLinear(x+1,y,s1)
                v[e] = 1
                e += 1  
            end 
            # vertical
            if y < s2
                i[e] = indexToLinear(x,y,s1)
                j[e] = indexToLinear(x,y+1,s1)
                v[e] = 1
                e += 1
            end 
        end
    end
    return IJV(n,nnz,i,j,v)
end

function ggrid2(s1::Integer,s2::Integer)
    A = sparse(ggrid2_ijv(s1,s2))
    return A + A'
end

function ggrid2_checkered_ijv(s1::Integer,s2::Integer,b1::Integer,b2::Integer,w::Real)
    l1 = div(s1,b1)
    l2 = div(s2,b2)

    n = s1*s2
    nnz = 2*n - s1 - s2
    
    i = Array{Int64,1}(undef,nnz)
    j = Array{Int64,1}(undef,nnz)
    v = Array{Float64,1}(undef,nnz)

    e = 1

    for x = 1:s1
        for y = 1:s2
            checkered = (mod(div(x-1,l1)+div(y-1,l2),2) == 0)
            # horizontal
            if x < s1
                i[e] = indexToLinear(x,y,s1)
                j[e] = indexToLinear(x+1,y,s1)
                v[e] = checkered ? w : 1
                e += 1
            end 
            # vertical
            if y < s2
                i[e] = indexToLinear(x,y,s1)
                j[e] = indexToLinear(x,y+1,s1)
                v[e] = checkered ? w : 1
                e += 1
            end 
        end
    end
    return IJV(n,nnz,i,j,v)
end

function ggrid2_checkered(s1::Integer,s2::Integer,b1::Integer,b2::Integer,w::Real)
    A = sparse(ggrid2_checkered_ijv(s1,s2,b1,b2,w))
    return A + A'
end

function ggrid2_checkered(s::Integer,b::Integer,w::Real)
    return ggrid2_checkered(s,s,b,b,w)
end

function plot_graph_weighted(gr,x,y;dots=true,setaxis=true,number=false)

    if isa(color, Vector) && length(color) == 3
        col = RGB(color...)
    else
        col = color
    end

    p = plot(;legend=false, axis=false, xticks=false, yticks=false)

    (ai,aj,av) = findnz(triu(gr))
    for i in 1:length(ai)
        s = [ai[i]; aj[i]]
        plot!(p, x[s], y[s],line_z=av[i], linecolor=:blues )
    end

    if dots
        scatter!(p,x, y , markercolor=col, markerstrokecolor=false)
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

    plot!(p; ylims = (miny - dely/20, maxy + dely/20))
    plot!(p; xlims = (minx - delx/20, maxx + delx/20))

    display(p)
    return p
end

function ggrid2coords(s1::Integer, s2::Integer)
    x = kron(ones(s2),collect(1:s1))
    y = kron(collect(1:s2),ones(s1))
    return x, y
end # grid2coords


function ggrid3_ijv(s1::Integer,s2::Integer,s3::Integer)
    n = s1*s2*s3
    nnz = 3*n - s1*s2 - s1*s3 - s2*s3
    
    i = Array{Int64,1}(undef,nnz)
    j = Array{Int64,1}(undef,nnz)
    v = Array{Float64,1}(undef,nnz)

    e = 1

    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                # eastward edge
                if x < s1
                    i[e] = indexToLinear(x,y,z,s1,s2)
                    j[e] = indexToLinear(x+1,y,z,s1,s2)
                    v[e] = 1
                    e += 1  
                end 
                # northward edge
                if y < s2
                    i[e] = indexToLinear(x,y,z,s1,s2)
                    j[e] = indexToLinear(x,y+1,z,s1,s2)
                    v[e] = 1
                    e += 1
                end 
                # upward edge
                if z < s3
                    i[e] = indexToLinear(x,y,z,s1,s2)
                    j[e] = indexToLinear(x,y,z+1,s1,s2)
                    v[e] = 1
                    e += 1
                end 
            end
        end
    end
    return IJV(n,nnz,i,j,v)
end

function ggrid3_checkered_ijv(s1::Integer,s2::Integer,s3::Integer,b1::Integer,b2::Integer,b3::Integer,w::Real)
    l1 = div(s1,b1)
    l2 = div(s2,b2)
    l3 = div(s3,b3)

    n = s1*s2*s3
    nnz = 3*n - s1*s2 - s1*s3 - s2*s3
    
    i = Array{Int64,1}(undef,nnz)
    j = Array{Int64,1}(undef,nnz)
    v = Array{Float64,1}(undef,nnz)

    e = 1

    for z = 1:s3
        for y = 1:s2
            for x = 1:s1
                checkered = (mod(div(x-1,l1)+div(y-1,l2)+div(z-1,l3),2) == 0)

                # eastward edge
                if x < s1
                    i[e] = indexToLinear(x,y,z,s1,s2)
                    j[e] = indexToLinear(x+1,y,z,s1,s2)
                    v[e] = checkered ? w : 1
                    e += 1  
                end 
                # northward edge
                if y < s2
                    i[e] = indexToLinear(x,y,z,s1,s2)
                    j[e] = indexToLinear(x,y+1,z,s1,s2)
                    v[e] = checkered ? w : 1
                    e += 1
                end 
                # upward edge
                if z < s3
                    i[e] = indexToLinear(x,y,z,s1,s2)
                    j[e] = indexToLinear(x,y,z+1,s1,s2)
                    v[e] = checkered ? w : 1
                    e += 1
                end 
            end
        end
    end
    return IJV(n,nnz,i,j,v)
end

function ggrid3(s1::Integer,s2::Integer,s3::Integer)
    A = sparse(ggrid3_ijv(s1,s2,s3))
    return A + A'
end

function ggrid3_checkered(s1::Integer,s2::Integer,s3::Integer,b1::Integer,b2::Integer,b3::Integer,w::Real)
    A = sparse(ggrid3_checkered_ijv(s1,s2,s3,b1,b2,b3,w))
    return A + A'
end

function ggrid3coords(s1::Integer, s2::Integer, s3::Integer)
    x = kron(ones(s3),kron(ones(s2),collect(1:s1)))
    y = kron(ones(s3),kron(collect(1:s2),ones(s1)))
    z = kron(collect(1:s3),kron(ones(s2),ones(s1)))
    return x, y, z
end # grid2coords

function plot_graph_weighted(gr,x,y,z)#;dots=true,setaxis=true)#,number=false)

    if isa(color, Vector) && length(color) == 3
        col = RGB(color...)
    else
        col = color
    end

    p = plot3d(;legend=false, axis=false, xticks=false, yticks=false, zticks=false)

    (ai,aj,av) = findnz(triu(gr))
    for i in 1:length(ai)
        s = [ai[i]; aj[i]]
        plot3d!(p, x[s], y[s], z[s], line_z=av[i], linecolor=:blues )
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

function getBoundary2(s1::Integer,s2::Integer)
    n = s1*s2
    bot = 1:s1
    top = n-s1+1:n
    left = s1+1:s1:n-s1
    right = 2*s1:s1:n-1
    bndry = [bot; top; left; right]
    return bndry
end

function getInterior2(s1::Integer,s2::Integer)
    n = s1*s2
    return setdiff(1:n,getBoundary2(s1,s2))
end

function getBoundary3(s1::Integer,s2::Integer,s3::Integer)
    n = s1*s2*s3
    zsidesize = s1*s2
    xsidesize = s2*s3
    ysidesize = s1*s3

    # alternatively, we could get z boundary as
    # zside1b = 1:zsidesize # bot 
    # zside2b = n-zsidesize+1:n # top
    
    zside1 = Array{Int64}(undef,zsidesize)
    for (i,(x,y)) in Base.Iterators.enumerate(Base.product(1:s1,1:s2))
        zside1[i] = indexToLinear(x,y,1,s1,s2)
    end
    
    zside2 = Array{Int64}(undef,zsidesize)
    for (i,(x,y)) in Base.Iterators.enumerate(Base.product(1:s1,1:s2))
        zside2[i] = indexToLinear(x,y,s3,s1,s2)
    end
    
    xside1 = Array{Int64}(undef,xsidesize)
    for (i,(y,z)) in Base.Iterators.enumerate(Base.product(1:s2,1:s3))
        xside1[i] = indexToLinear(1,y,z,s1,s2)
    end
    
    xside2 = Array{Int64}(undef,xsidesize)
    for (i,(y,z)) in Base.Iterators.enumerate(Base.product(1:s2,1:s3))
        xside2[i] = indexToLinear(s1,y,z,s1,s2)
    end
    
    yside1 = Array{Int64}(undef,ysidesize)
    for (i,(x,z)) in Base.Iterators.enumerate(Base.product(1:s1,1:s3))
        yside1[i] = indexToLinear(x,1,z,s1,s2)
    end
    
    yside2 = Array{Int64}(undef,ysidesize)
    for (i,(x,z)) in Base.Iterators.enumerate(Base.product(1:s1,1:s3))
        yside2[i] = indexToLinear(x,s2,z,s1,s2)
    end
    
    return [zside1; zside2; xside1; xside2; yside1; yside2]
end

function getInterior3(s1::Integer,s2::Integer,s3::Integer)
    n = s1*s2*s3
    return setdiff(1:n,getBoundary3(s1,s2,s3))
end

# function zeroBoundaryLap2(s1,)