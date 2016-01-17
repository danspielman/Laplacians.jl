% add current directory to the path
addpath(pwd())

% if Julia is not on the search path, try to find it.

status = system('julia --version');

if (status ~= 0)
    if (isunix)
        
        % if we can find /Applications, then find most recent
        % version of Julia
        
         strs = dir('/Applications/Julia*');
         
         vern = zeros(length(strs),1);
         
         for i = 1:length(strs)
           str = strs(i).name;
           ind = (str >= '0') .* (str <= '9');
           vern(i) = str2num(['1', str(find(ind))]);
         end
         [~,i] = max(vern);
         
         p = ['/Applications/', strs(i).name, '/Contents/', ...
                             'Resources/julia/bin/'];
              
         setenv('PATH',[':', p]);
        
    else
        error('Cannot find Julia in the shell path')
    end
end

status = system('julia --version');
if (status == 0)
    fprintf('Found Julia and put it on the system path.\n')
end

    