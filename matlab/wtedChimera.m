function a = wtedChimera(n, i)
% function a = wtedChimera(n, i)
%
% produce the ith weighted chimera graph on n vertices from Lapsolver.jl,
% using a call to julia

nm = [tempname, '.txt'];

cmd = ['julia '];

cmd = [cmd, ' -e ''using Laplacians; writeIJV("', nm, '",' ...
                    'wtedChimera('];
       
cmd = [cmd, num2str(n), ',', num2str(i) '))'' '];


system(cmd);

a = readIJV(nm);
delete(nm);
