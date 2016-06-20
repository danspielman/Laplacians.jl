function x = lapCmgSolver(la,b,opts)
% function x = lapCmgSolver(la,b,opts)
%
% wrapping of lapWrapSolver around cmgSolver
%

default('b',[]);

default('opts','tol',1e-6);
default('opts','maxit',100);

if (isempty(b)),
    f = lapWrapSolver('cmgSolver',la,[],opts);
    x = f;
    return
else
  x = lapWrapSolver('cmgSolver',la,b,opts);
end

    
