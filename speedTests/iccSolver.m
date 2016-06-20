function [x] = iccSolver(la,b,opts)
% function [x] = iccSolver(la,b,opts)
% function [f] = iccSolver(la,[],opts)
%
% this calls pcg with incomplete cholesky preconditioner,
% opts is passed to ichol, defaults are:
%   'nofill' and michol 'off'
%   tol 1e-6
%   maxit 100
%   if opts.L is supplied, it skips the call to ichol
%
% puts matrix in rcm order
%
%
    
%
% to do: allow other orders


default('b',[]);
default('opts','type','nofill');
default('opts','tol',1e-6);
default('opts','maxit',10000);
default('opts','L',[]);

%opts.michol = 'on';

p = symrcm(la);
laperm = la(p,p);

icholOpts.type = opts.type;

if (isempty(opts.L))
  L2 = ichol(laperm,icholOpts);
else
  L2 = opts.L;
end

if isempty(b)
    f = @(b)(internal(laperm,p,L2,b,opts));
    x = f;
else
    x = internal(laperm,p,L2,b,opts);
end

end % main function

function x = internal(laperm,p,L2,b,opts)

  bperm = b(p);
  [xperm,flag,relres,iter] = pcg(laperm,bperm,opts.tol,opts.maxit,L2,L2');
  x(p) = xperm;
  x = x(:);

end
  


