immutable elimTreeNode
    nodeid::Int64
    parent::Int64
    wtDeg::Float64
end

function makeElimList(t)
    tr = matToTree(t)
    n = size(tr.children,1)  
    
    elims = Array{elimTreeNode}(0)
    for vi in n:-1:2
        v = tr.children[vi]
        push!(elims,elimTreeNode(v,tr.parent[v],tr.weights[vi]))
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