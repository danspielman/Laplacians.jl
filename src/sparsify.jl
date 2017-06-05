#=

Code for Sparsification.

Right now, just implements Spielman-Srivastava

=#

function sparsify(a; ep=0.3, matrixConcConst=4.0, JLfac=16)
  f = approxCholLap(a,tol=1e-2);

  n = size(a,1)
  @show k = round(Int, JLfac*log(n)) # number of dims for JL

  U = wtedEdgeVertexMat(a)
  m = size(U,1)
  R = randn(m,k)
  UR = U'*R;

  V = zeros(n,k)
  for i in 1:k
    V[:,i] = f(UR[:,i])
  end

  (ai,aj,av) = findnz(triu(a))
  prs = zeros(av)
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
