function writeIJV(filename, a)
% function writeIJV(filename, a)
%
% write upper-triangular portion of a as a comma-separated IJV

[ai,aj,av] = find(triu(a));
dlmwrite(filename,[ai,aj,av],'precision',9);
