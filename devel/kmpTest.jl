
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



function testSampler(a; t=akpw(a), frac1=1/5)


    n = size(a,1);

    rest = a-t;

    st = compStretches(t,rest);
    aveStretch = sum(st)/nnz(rest)
    @show aveStretch
    

    targetStretch = 1/(2*log(n)/log(2))

    fac = aveStretch/targetStretch
    heavy = rest+fac*t;

    (ai,aj,av) = findnz(triu(rest))
    (si,sj,sv) = findnz(triu(st))
    sv = sv ./ fac

    ijvs = IJVS(ai,aj,av,sv)

    ijvs1 = stretchSample(ijvs,targetStretch,frac1)
    ijvs2 = stretchSample(ijvs1,targetStretch,frac1)
    
    samp1 = sparse(ijvs1.i,ijvs1.j,ijvs1.v,n,n)
    samp1 = samp1 + samp1';
    
    samp2 = sparse(ijvs2.i,ijvs2.j,ijvs2.v,n,n)
    samp2 = samp2 + samp2';
    
    add = speye(n)/10^6;
    e1 = maximum(eigs(lap(heavy)+add,lap(samp1+t*fac)+add, nev=1)[1]);
    e2 = maximum(eigs(lap(samp1+t*fac)+add,lap(samp2+t*fac)+add, nev=1)[1]);

    return [e1,e2]
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

# makes the spine heavy graph, to test out the preconditioner
function makeHeavy(a; t=akpw(a),  params::KMPparams=defaultKMPparams)
   n = size(a,1);

    ord = Laplacians.dfsOrder(t)

    aord = a[ord,ord]
    tord = t[ord,ord]

    rest = aord-tord;

    st = compStretches(tord,rest);
    aveStretch = sum(st)/nnz(rest)

    targetStretch = 1/(params.treeScale*log(n)/log(2))

    fac = aveStretch/targetStretch
    tree = fac*t;

    heavy = rest + tree;

    return heavy
end



function vecstats(s)
    println("length : ", size(s,1), ", min : ", minimum(s), ", mean : ", mean(s), ", max : ", maximum(s), ", sum : ", sum(s))
end


"""Do `numruns` tests on chimeras of size `n`, going through each param choice
in `pList`.  Also use augTreeSolver, for Laplacians, for comparison (last).
Return the times of all runs. For example:

~~~
pList = [KMPparams(1/6^2,6,1/4,600,:rand)]
push!(pList,KMPparams(1/4^2,4,1/4,600,:rand))
out = manyRuns(10000,10,pList)
~~~
"""
function manyRuns(n,numruns,pList)

    out = zeros(numruns, 1+length(pList))
    tot = zeros(1+length(pList))
    for it in 1:numruns
        println("chimera(", n, ", ", it, ")" )
        a = chimera(n,it)
        b = randn(n)
        b = b - mean(b)
        
        for i in 1:length(pList)
            fsub =  KMPLapSolver(a, tol = 0.01, maxits = 10000, params=pList[i], verbose=false)
            tic()
            xh = fsub(b)
            ti = toq()
            tot[i] += ti
            out[it,i] = ti
        end

        fsub = lapWrapSolver(augTreeSolver,lap(a),tol=0.01,maxits=10000)
        tic()
        xh = fsub(b)
        ti = toq()
        tot[1+length(pList)] += ti
        out[it,1+length(pList)] = ti


        for i in 1:length(pList)
            println(string(pList[i]), " : ", tot[i])
        end
        println("augTree : ", tot[1+length(pList)])

    end

            

    return out
end

function manyRunsW(n,numruns,pList)

    out = zeros(numruns, 1+length(pList))
    tot = zeros(1+length(pList))
    for it in 1:numruns
        println("wtedChimera(", n, ", ", it, ")" )
        a = wtedChimera(n,it)
        b = randn(n)
        b = b - mean(b)
        
        for i in 1:length(pList)
            fsub =  KMPLapSolver(a, tol = 0.01, maxits = 10000, params=pList[i], verbose=false)
            tic()
            xh = fsub(b)
            ti = toq()
            tot[i] += ti
            out[it,i] = ti
        end

        fsub = lapWrapSolver(augTreeSolver,lap(a),tol=0.01,maxits=10000)
        tic()
        xh = fsub(b)
        ti = toq()
        tot[1+length(pList)] += ti
        out[it,1+length(pList)] = ti


        for i in 1:length(pList)
            println(string(pList[i]), " : ", tot[i])
        end
        println("augTree : ", tot[1+length(pList)])

    end

            

    return out
end
