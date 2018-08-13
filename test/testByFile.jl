# Tests for functions listed by file
# 

# solverInterface.jl

a = randn(10,10)
a = Laplacians.forceLap(a)

Laplacians.forceLap(a)
Laplacians.forceLap(abs.(a))


a = wted_chimera(100,101);
la = lap(a);
la[1,1] += 1;
f = Laplacians.sddmWrapLap(approxCholLap, la, verbose=true);
b = randn(100);
x = f(b);

a = grid2(20)
a2 = subsampleEdges(grid2(20),0.45)
f = approxCholLap(a2)
its = [0]
f = approxCholLap(a2, pcgIts=its, verbose=true)

# conditionNumber

support(a,a2)
support(a2,a)

conditionNumber(a, akpw(a), verbose=true)

f = approxCholLap(a, maxits=1)
conditionNumber(a, f, verbose=true)

println("End of testByFile")
