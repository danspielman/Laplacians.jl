function lexCoords(a,x,y; tol=.0001, maxits=1000)

    n = a.n
    
    maxx = zeros(n)
    minx = zeros(n)
    maxy = zeros(n)
    miny = zeros(n)
    
    iters = 0
    while iters < maxits
    
        x = x - mean(x)
        y = y - mean(y)
        
        for v in 1:n
            maxx[v] = maximum(x[a.rowval[a.colptr[v]:(a.colptr[v+1]-1)]])
            minx[v] = minimum(x[a.rowval[a.colptr[v]:(a.colptr[v+1]-1)]])
            maxy[v] = maximum(y[a.rowval[a.colptr[v]:(a.colptr[v+1]-1)]])
            miny[v] = minimum(y[a.rowval[a.colptr[v]:(a.colptr[v+1]-1)]])
        end
        
        prevx = x
        
        for v in 1:n
            x[v] = (2*x[v] + maxx[v]+minx[v])/4
            y[v] = (2*y[v] + maxy[v]+miny[v])/4
        end
        
        xy = [x y];
        fix = inv(sqrtm(xy'*xy))
        xy = xy*fix
        x = xy[:,1]
        y = xy[:,2]
        
        if (norm(prevx-x)<norm(x)*tol)
            break
        end
        
        iters = iters+1

    end
    
    return x,y
end

function lexStar(v, w)
    if v[1] < v[2]
        low = 1
        high = 2
    else
        low = 2
        high = 1
    end
    del = w[low]*w[high]/(w[low]+w[high])*(v[high]-v[low])
    x = (w[low]*v[low] + w[high]*v[high])/(w[low]+w[high])

    for i in 3:length(v)
        deli = w[i]*abs(v[i] - x)
        if (deli > del)
            
            # pick the best other end
            other = 1
            del = w[other]*w[i]/(w[other]+w[i])*abs(v[other]-v[i])
            for j in 2:(i-1)
                delj = w[j]*w[i]/(w[j]+w[i])*abs(v[j]-v[i])
                if delj > del
                    other = j
                    del = delj
                end            
            end
            if v[other] < v[i]
                low = other
                high = i
            else
                low = i
                high = other
            end
            x = (w[low]*v[low] + w[high]*v[high])/(w[low]+w[high])
        end
    end
    return x
end

