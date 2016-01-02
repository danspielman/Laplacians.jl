immutable elimLeafNode
    nodeid::Int64
    parent::Int64
    wtDeg::Float64
end

immutable elimDeg2Node
    nodeid::Int64
    nbr1::Int64
    nbr2::Int64
    wt1::Float64
    wt2::Float64
end


function makeElimList(t)
    tr = matToTree(t)
    n = size(tr.children,1)  
    
    elims = Array{elimLeafNode}(0)
    for vi in n:-1:2
        v = tr.children[vi]
        push!(elims,elimLeafNode(v,tr.parent[v],tr.weights[vi]))
    end
    
    return elims
end

function solveTree(elims,b)
    cumb = copy(b)
    n = size(b,1)
    for i in 1:(n-1)
        cumb[elims[i].parent] += cumb[elims[i].nodeid]
    end
    x = zeros(Float64,n)
    for i in (n-1):-1:1
        node = elims[i].nodeid
        x[node] = x[elims[i].parent] + cumb[node]/elims[i].wtDeg
    end
    return x
    
end

# Must be in DFS order
# marked is 1 if flagged for possible elimination,
# and set to 2 if we do eliminate it
function elimDeg1(t, marked)
    n = t.n

    deg = Array{Int64}(n)
    for v in 1:n
        deg[v] = t.colptr[v+1] - t.colptr[v]
    end

    elims1 = Array{elimLeafNode}(0)

    for v in n:-1:2

        if (deg[v] == 1 && marked[v] == 1)
            parent = t.rowval[t.colptr[v]];
            wt = t.nzval[t.colptr[v]];
            push!(elims1,elimLeafNode(v,parent,wt))

            deg[parent] = deg[parent] - 1
            marked[v] = 2
        end
    end
    return elims1
end


# Must be in DFS order
# marked is 1 if flagged for possible elimination,
# and set to 2 if we do eliminate it
function elimDeg12(t, marked)

    # make sure root is not marked
    marked[1] = 0

    deg = Array{Int64}(n)
    for v in 1:n
        deg[v] = t.colptr[v+1] - t.colptr[v]
    end

    elims1 = Array{elimLeafNode}(0)

    for v in n:-1:2

        if (deg[v] == 1 && marked[v] == 1)
            parent = t.rowval[t.colptr[v]];
            wt = t.nzval[t.colptr[v]];
            push!(elims1,elimLeafNode(v,parent,wt))

            deg[parent] = deg[parent] - 1
            marked[v] = 2
            deg[v] = 0
        end
    end

    elims2 = Array{elimDeg2Node}(0)

    subt = triu(t)

    for v in n:-1:2

        if (deg[v] == 2 && marked[v] == 1)

            parent = t.rowval[t.colptr[v]];

            # to ident the child, enumerate to find one uneliminated, which we check by marked
            kidind = t.colptr[v]+1
            kid = t.rowval[kidind]
            while deg[kid] == 0
                kidind = kidind+1
                kid = t.rowval[kidind]

                if kidind >= t.colptr[v+1]
                    error("went of the end of the kid list without finding node to elim from")
                end
            end


            wt1 = t.nzval[t.colptr[v]];
            wt2 = t.nzval[kidind]

            push!(elims2,elimDeg2Node(v,parent,kid,wt1,wt2))
            marked[v] = 2

            newwt = 1/(1/wt1 + 1/wt2)

            # now that we've found the kid, go up the chain until done
            while (deg[parent] == 2 && marked[parent] == 1)
                v = parent
                parent = t.rowval[t.colptr[v]];
                wt1 = t.nzval[t.colptr[v]];
                wt2 = newwt

                push!(elims2,elimDeg2Node(v,parent,kid,wt1,wt2))
                marked[v] = 2
                
                newwt = 1/(1/wt1 + 1/wt2)
            end

            # now, hack the tree to adjust parent and wt of node kid
            subt.rowval[subt.colptr[kid]] = parent
            subt.nzval[subt.colptr[kid]] = newwt

        end
    end

    subt = subt + subt'

    ind = find(marked.<2)
    subt = subt[ind,ind]
    
    return elims1, elims2, ind, subt
end


function recTest1(t,b, marked)

    n = length(b)

    elims1 = elimDeg1(t, marked)

    cumb = copy(b)

    for i in 1:length(elims1)
        cumb[elims1[i].parent] += cumb[elims1[i].nodeid]
    end

    ind = find(marked.<2)
    cumbs = cumb[ind]

    @show sum(cumbs)

    ts = t[ind,ind];
    lts = lap(ts)

    xs = pinv(full(lts))*cumbs
    
    @show norm(lts*xs-cumbs)

    x = zeros(n)
    x[ind] = xs
    for i in length(elims1):-1:1
        node = elims1[i].nodeid
        x[node] = x[elims1[i].parent] + cumb[node]/elims1[i].wtDeg
    end

    return x
    
end


function recTest12(t,b, marked)

    n = length(b)

    elims1, elims2, ind, subt = elimDeg12(t, marked)

    @show length(elims1)
    #@show elims2

    cumb = copy(b)

    for i in 1:length(elims1)
        cumb[elims1[i].parent] += cumb[elims1[i].nodeid]
    end

    for i in 1:length(elims2)
        wtsum = elims2[i].wt1 + elims2[i].wt2
        cumb[elims2[i].nbr1] += cumb[elims2[i].nodeid]*elims2[i].wt1 / wtsum
        cumb[elims2[i].nbr2] += cumb[elims2[i].nodeid]*elims2[i].wt2 / wtsum
    end

    
    cumbs = cumb[ind]

    @show sum(cumbs)

    # this is not the right elim for deg2
    
    lts = lap(subt)

    xs = pinv(full(lts))*cumbs
    
    @show norm(lts*xs-cumbs)

    x = zeros(n)
    x[ind] = xs

    for i in length(elims2):-1:1
        node = elims2[i].nodeid
        wtsum = elims2[i].wt1 + elims2[i].wt2

        x[node] = (elims2[i].wt1*x[elims2[i].nbr1] + elims2[i].wt2*x[elims2[i].nbr2] + cumb[node])/wtsum
    end

    
    for i in length(elims1):-1:1
        node = elims1[i].nodeid
        x[node] = x[elims1[i].parent] + cumb[node]/elims1[i].wtDeg
    end

    return x
    
end



function forwardSolve(b, elims1, elims2)

    y = copy(b)

    for i in 1:length(elims1)
        y[elims1[i].parent] += y[elims1[i].nodeid]
    end

    for i in 1:length(elims2)
        wtsum = elims2[i].wt1 + elims2[i].wt2
        y[elims2[i].nbr1] += y[elims2[i].nodeid]*elims2[i].wt1 / wtsum
        y[elims2[i].nbr2] += y[elims2[i].nodeid]*elims2[i].wt2 / wtsum
    end

    return y
    
end

function backSolve(x, y, elims1, elims2)
    
    for i in length(elims2):-1:1
        node = elims2[i].nodeid
        wtsum = elims2[i].wt1 + elims2[i].wt2

        x[node] = (elims2[i].wt1*x[elims2[i].nbr1] + elims2[i].wt2*x[elims2[i].nbr2] + y[node])/wtsum
    end

    
    for i in length(elims1):-1:1
        node = elims1[i].nodeid
        x[node] = x[elims1[i].parent] + y[node]/elims1[i].wtDeg
    end

end


function recTest12x(t,b, marked)

    n = length(b)

    elims1, elims2, ind, subt = elimDeg12(t, marked)

    lts = lap(subt)
    lts[1,1] = lts[1,1] + 1
    
    F = factorize(lts)


    y = forwardSolve(b, elims1, elims2)

    ys = y[ind]


    xs = F \ ys
    
    x = zeros(n)
    x[ind] = xs

    backSolve(x, y, elims1, elims2)
    
    return x
    
end






function testKMP(a; t=akpw(a), frac1=1/5, frac2=1/20)
    st = compStretches(t,a);
    aveStretch = sum(st)/nnz(a)
    @show aveStretch
    
    rest = a-t;
    fac = aveStretch*log(size(a,1))*3
    heavy = rest+fac*t;
    n = size(a,1);
    
    strest = compStretches(t,rest)
    probs = triu(strest);
    probs = probs * (n*frac1)/ sum(probs);
    
    (pi,pj,pv) = findnz(probs)
    select = rand(size(pv)) .< pv;
    (ai,aj,av) = findnz(triu(rest))
    
    samp1 = sparse(ai,aj,av.*select./pv,n,n)
    samp1 = samp1 + samp1';
    
    st1 = compStretches(t*fac, samp1);
    (s1i,s1j,s1v) = findnz(triu(samp1))
    
    nits = 20
    dat = zeros(nits)
    for it in 1:nits
        pr = (frac2/frac1)
        select1 = rand(size(s1v)) .< pr
        samp2 = sparse(s1i,s1j,s1v.*select1,n,n)
        samp2 = samp2 + samp2';
        add = speye(n)/10^6;

        e = eigs(lap(samp1+t*fac)+add,lap(samp2+t*fac)+add, nev=4);

        dat[it] = maximum(e[1])
    end
    [minimum(dat) , mean(dat), median(dat), maximum(dat)]
end

function testAtSize(n, results; frac1=1/3, frac2 = 1/15)
    mxval = 0
    mxi = 0

    i = 1
    while i < 10^6
        a = wtedChimera(n,i)
        x = testKMP(a, frac1 = frac1, frac2 = frac2)
        if maximum(x) > mxval
            mxi = i
            mxval = maximum(x)
        end

        println("on iter " , i , ". Max was ", mxi, " : ", mxval )

	push!(results,maximum(x))

        i = i + 1
    end

    
end