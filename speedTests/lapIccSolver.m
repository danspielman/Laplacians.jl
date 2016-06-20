function x = lapIccSolver(la,b,opts)
% function x = lapIccSolver(la,b,opts)
%
% wrapping of lapWrapSolver around iccSolver
%

default('opts','type','nofill');
default('b',[]);

if (isempty(b)),
    f = lapWrapSolver('iccSolver',la,[],opts);
    x = f;
    return
else
  x = lapWrapSolver('iccSolver',la,b,opts);
end

    
