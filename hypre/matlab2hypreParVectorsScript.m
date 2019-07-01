% should be called with 
% EXAMPLE
%  /Applications/MATLAB_R2018b.app/bin/matlab -nojvm < matlab2hypreParVectorsScript" 

display(pwd)

load 'julia2matlab2hypre_vector.mat'
% hardcoded vector file name, 
% must load
%   vector with name b
%   string with name output_filename
%   integer with name num_procs

filename = output_filename;

precision = '16.15e';
prt_format = strcat('%', precision);

[n, m] = size(b);

% Generate partitioning across processors.
% See Hypre getpart.c.
part_size = floor(n/num_procs);
rest = mod(n, num_procs);
part = [0 (rest + part_size):part_size:n];

% generate Hypre input files
Y = [n part(1:num_procs)]'; %comment here just to kill bug in vscode syntax highlighting by closing transpose that it thinks starts a string'

for i = 1:num_procs
    nrows = part(i+1) - part(i);
    filename2 = [filename, '.', num2str(i-1)];
    fprintf('Generating file: %s\n', filename2);
    X = b((part(i) + 1):part(i+1));

    dlmwrite(filename2, nrows, 'precision', '%d');
    dlmwrite(filename2, X, '-append', 'precision', prt_format);

    % writing INFO file
    filename2 = [filename, '.', 'INFO', '.',...
        num2str(i-1)];
    fprintf('Generating INFO file: %s\n', filename2);
    dlmwrite(filename2, Y, 'precision', '%d');
end

exit