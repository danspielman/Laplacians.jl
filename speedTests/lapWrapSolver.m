function x = lapWrapSolver(solvername,la,b,opts)
% function x = lapWrapSolver(solvername,la,b,opts)
% function f = lapWrapSolver(solvername,la,[],opts)
% function f = lapWrapSolver(solvername,la)
% function f = lapWrapSolver(solvername)
%
% takes a solver for pos def systems,
% and uses it to construct a solver for a laplacian system
% by fixing last entry
%
% if only solvername is supplied, this creates a wrapped version of
% that solver.  But, it can only be applied to a matrix.
%
% if the solver takes options, it can supply those too
%
% Examples:
%  a = grid2(10); 
%  la = lap(a);
%  b = randn(length(la),1); b = b - mean(b);
%  x = lapWrapSolver('iccSolver',la,b); norm(la*x-b)
%  f = lapWrapSolver('iccSolver',la); norm(la*f(b)-b)
%  g = lapWrapSolver('iccSolver'); f = g(la); norm(la*f(b)-b)  
%
% Copyright Daniel Spielman, 2013, Yale University.
% This file of the LapSolve package.
% See source file for license information

% This file is part of LapSolve.
%
% LapSolve is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% LapSolve is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with LapSolve.  If not, see <http://www.gnu.org/licenses/>.

default('b',[]);
default('opts',[]);

if (nargin==1),
    f = @(la)(lapWrapSolver(solvername,la,b,opts));
    x = f;
    return
end

    

n = length(la);
lasub = la(1:(n-1),1:(n-1));
          
if isempty(b)
  
  str = ['@(mat,b,opts)(', solvername, '(mat,b,opts))'];
  g = eval(str);
  g2 = g(lasub,[],opts);
  

  f = @(b)(internal(g2,n,b));
  x = f;

else
  
  str = ['@(mat,b,opts)(', solvername, '(mat,b,opts))'];
  g = eval(str);
  g2 = g(lasub,[],opts);
  
  x = internal(g2,n,b);

end

end % main function
          
          
function x = internal(g2,n,b)

    b = b - mean(b);
    x = g2(b(1:(n-1)));
    x(n) = 0;
    x = x - mean(x);
    x = x(:);
    
end

