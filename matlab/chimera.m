function a = chimera(n, i)
% function a = chimera(n, i)
%
% produce the ith chimera graph on n vertices from Lapsolver.jl,
% using a call to julia

nm = [tempname, '.txt'];

cmd = ['julia '];

cmd = [cmd, ' -e ''using Laplacians; writeIJV("', nm, '",' ...
                    'chimera('];
       
cmd = [cmd, num2str(n), ',', num2str(i) '))'' '];

system(cmd)

a = readIJV(nm);
delete(nm);
