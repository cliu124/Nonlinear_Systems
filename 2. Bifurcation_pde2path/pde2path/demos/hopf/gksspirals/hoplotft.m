function hoplotft(dir,pt,wnr,cnr,varargin)
% HOPLOTF: call hoplot for point read from dir/pt
%
%  hoplotf(dir,pt,wnr,cnr,varargin)
%
% varargin passed through to hopot
try p=loadp(dir,pt); 
catch fprintf(['point ' dir '/' pt ' does not exist.\n']); 
    npt=ptselect(dir); p=loadp(dir,npt); pt=npt; 
end
fprintf('lam=%g, T=%g\n',getlam(p),p.hopf.T); 
if nargin>4; hoplott(p,wnr,cnr,varargin{1}); 
else try; hoplott(p,wnr,cnr,p.hopf.aux); catch; hoplot(p,wnr,cnr);   end 
end
