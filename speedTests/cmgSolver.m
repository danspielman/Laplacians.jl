function [x] = cmgSolver(la,b,opts)
% function [x] = cmgSolver(la,b,opts)
% function [f] = cmgSolver(la,[],opts)
%
% this calls Yiannis Koutis's cmg routines.
% they must be in the default path to work
%
% If la is singular, wrap with lapWrapSolver
%
% opts.tol (default 1e-6)
% opts.maxit (default 100)
%


default('b',[]);

default('opts','tol',1e-6);
default('opts','maxit',100);

cmgOpts.display = 0;
pfun = cmg_sdd(la,cmgOpts);

if isempty(b)
    f = @(b)(internal(la,pfun,b,opts));
    x = f;
else
    x = internal(la,pfun,b,opts);
end

end % main function


function x = internal(la,pfun,b,opts)

  [x,flag] = pcg(la, b, opts.tol, opts.maxit, pfun);
  %  [x] = pcg(la, b, opts.tol, opts.maxit, pfun);

end
