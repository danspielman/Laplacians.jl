#=

Code for Sparsification.

Right now, just implements Spielman-Srivastava

=#

"""
    as = sparsify(a; ep=0.5)

Apply Spielman-Srivastava sparsification: sampling by effective resistances.
`ep` should be less than 1.
"""
function sparsify(a; ep=0.3, matrixConcConst=4.0, JLfac=4.0)

  if ep > 1
    @warn "Calling sparsify with ep > 1 can produce a disconnected graph."
  end

  f = approxCholLap(a,tol=1e-2);

  n = size(a,1)
  k = round(Int, JLfac*log(n)) # number of dims for JL

  U = wtedEdgeVertexMat(a)
  m = size(U,1)
  R = randn(m,k)
  UR = U'*R;

  V = zeros(n,k)
  for i in 1:k
    V[:,i] = f(UR[:,i])
  end

  (ai,aj,av) = findnz(triu(a))
  prs = zeros(size(av))
  for h in 1:length(av)
      i = ai[h]
      j = aj[h]
      prs[h] = min(1,av[h]* (norm(V[i,:]-V[j,:])^2/k) * matrixConcConst *log(n)/ep^2)
  end

  ind = rand(size(prs)) .< prs

  as = sparse(ai[ind],aj[ind],av[ind]./prs[ind],n,n)
  as = as + as'

  return as

end
