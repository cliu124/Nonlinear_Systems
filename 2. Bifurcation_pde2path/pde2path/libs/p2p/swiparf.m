function p=swiparf(dir,ptnr,ndir,par)
% SWIPARF: switch parameter, solution from file
%
%  p=swiparf(dir,fname,ndir,par)
%
% use p-struct in dir/fname and pass ndir, par to swipar
%
% See also swipar, loadp
try p=loadp(dir,ptnr); 
catch nptnr=char(['pt' mat2str(max(getlabs(dir)))]); 
    fprintf(['point ' dir '/' ptnr ' does not exist. Please check what points are saved in ' dir '.\n']); 
    fprintf(['Instead using ' dir '/' nptnr '.\n']);
    p=loadp(dir,nptnr); 
end 
p=swipar(p,par,ndir);
end
