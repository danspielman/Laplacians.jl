a = wtedChimera(2000)
m = randperm(nnz(a))
list = m[1:1000]
list = sort(list);

b = submatrixCSC(a,list)
(ai,aj,av) = findnz(a)
b2 = sparse(ai[list],aj[list],av[list],a.m,a.n)

@test sum(abs(b - b2)) == 0


perm = randperm(2000)
ap = symPermuteCSC(a, perm)
ap2 = a[perm,perm]
@test sum(abs(ap-ap2)) == 0


