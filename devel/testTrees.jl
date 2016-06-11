function testTrees(n,nruns,fn)

st0 = Array(Float64,0)
st2 = Array(Float64,0)
st4 = Array(Float64,0)
st5 = Array(Float64,0)
stp = Array(Float64,0)
stk = Array(Float64,0)
stold = Array(Float64,0)

t0 = Array(Float64,0)
t2 = Array(Float64,0)
t4 = Array(Float64,0)
t5 = Array(Float64,0)
tp = Array(Float64,0)
tk = Array(Float64,0)
told = Array(Float64,0)

h = open(fn,"w")

for i in 1:nruns
    a = wtedChimera(n,i)
    f(t) = sum(compStretches(t,a))/nnz(a)
    
    try
    
    tic()
    tr = akpwish(a,ver=0)
    at0 = toq()
    ast0 = f(tr)

    tic()
    tr = akpwish(a,ver=2)
    at2 = toq()
    ast2 = f(tr)

    tic()
    tr = akpwish(a,ver=4)
    at4 = toq()
    ast4 = f(tr)

        #=
    tic()
    tr = akpwish(a,ver=5)
    at5 = toq()
    ast5 = f(tr)
        =#
        
    tic()
    tr = randishKruskal(a)
    atk = toq()
    astk = f(tr)
    
    tic()
    tr = randishPrim(a)
    atp = toq()
    astp = f(tr)
    
    tic()
    tr = akpw(a)
    atold = toq()
    astold = f(tr)
    
    print(h, n, " ", i, " : ")
    
    print(h, ast0, " ")
    push!(t0,at0)
    push!(st0,ast0)

    print(h, ast2, " ")
    push!(t2,at2)
    push!(st2,ast2)
    
    print(h, ast4, " ")
    push!(t4,at4)
    push!(st4,ast4)
    
        #=
    print(h, ast5, " ")
    push!(t5,at5)
    push!(st5,ast5)
        =#
        
    print(h, astk, " ")
    push!(tk,atk)
    push!(stk,astk)
    
    print(h, astp, " ")
    push!(tp,atp)
    push!(stp,astp)
    
    print(h, astold, " ")
    push!(told,atold)
    push!(stold,astold)
    
    print(h, "\n")
        
    catch
        println("bad on wtedChimera ", n, " ", i)
    end


end
       


mist = minimum([st0 st2 st4 stk stp stold],2)
mint = minimum([t0 t2 t4  tk tp told],2)


f = 1.3

stats(v) = [mean(v./mist) median(v./mist) maximum(v./mist) mean(v .== mist) mean(v .> f*mist)]
statt(v) = [mean(v./mint) median(v./mint) maximum(v./mint) mean(v .== mint) mean(v .> f*mint)]

println(h, "st0 : ", stats(st0))
println(h, "st2 : ", stats(st2))
println(h, "st4 : ", stats(st4))
#println(h, "st5 : ", stats(st5))
println(h, "stp : ", stats(stp))
println(h, "stk : ", stats(stk))
println(h, "stold : ", stats(stold))

println(h, "t0 : ", statt(t0))
println(h, "t2 : ", statt(t2))
println(h, "t4 : ", statt(t4))
#println(h, "t5 : ", statt(t5))
println(h, "tp : ", statt(tp))
println(h, "tk : ", statt(tk))
println(h, "told : ", statt(told))

    close(h)

end
