# totTime should be in hours
function compareSolversOut(n, solvers, tol, totTime)
    
    # first, force a compile
    tab = compareSolversRuns(1000, solvers, tol, 2)

    tab = compareSolversTime(n, solvers, tol, totTime)
    nruns = size(tab,1) - 1

    fn = string("compSolvers", n , "x", nruns, ".csv")
    h = open(fn,"w")
    write(h, chomp(readall(`hostname`)))
    write(h, "  ")
    write(h, string(now()))
    write(h,"\n")

    appendcsv(fn,tab)

end


compareSolversRuns(n, solvers, tol, nruns, maxtime=Inf) = compareSolvers(n, solvers, tol=tol, nruns=nruns, maxtime=maxtime)
compareSolversTime(n, solvers, tol, totTime) = compareSolvers(n, solvers, tol=tol, totTime=totTime)

# totTime comes in hours, convert to seconds
function compareSolvers(n, solvers; tol=1e-3, nruns=10^8, totTime=Inf, maxtime=totTime*60*60/10)

    totTime = totTime*60*60

    numSolvers = length(solvers)
    bt = []
    st = []
    acc = []
    name = Array{ASCIIString,1}(0)
    for i in 1:numSolvers
        push!(bt, Array{Float64,1}(0))
        push!(st, Array{Float64,1}(0))
        push!(acc, Array{Float64,1}(0))
        push!(name, split(string(solvers[i]),'.')[2])
    end
    nedges = Array(Int64,0)
    nlist = Array(Int64,0)
    ilist = Array(Int64,0)


    t1 = time()
    for run in 1:nruns
        if time() > t1 + totTime
            break
        end

        a = wtedChimera(n, run)
        la = lap(a)

        b = randn(n)
        b = b - mean(b)

        push!(nlist, n)
        push!(ilist, run)
        push!(nedges, nnz(a))

        for i in 1:numSolvers

            tic()
            f = solvers[i](a, tol=tol, maxtime=maxtime)
            push!(bt[i], toq());

            tic()
            x = f(b);
            push!(st[i], toq());

            push!(acc[i], norm(la * x - b) / norm(b))

        end

        print(".")

    end

    tab = ["nlist"; nlist]
    tab = [tab ["ilist"; ilist]]
    tab = [tab ["nedges"; nedges]]

    for i in 1:numSolvers
        tab = [tab [name[i] * "Build"; bt[i]]]
        tab = [tab [name[i] * "SolveTime"; st[i]]]
        tab = [tab [name[i] * "RelativeError"; acc[i]]]
    end

    return tab
end

function appendcsv(filename, data)

    (n1,n2) = size(data)
    
    fh = open(filename,"a")
    for i in 1:n1
        for j in 1:n2
            write(fh, string(data[i,j]))
            if j < n2
                write(fh, ", ")
            end
        end
        write(fh, "\n")
    end

    close(fh)
end