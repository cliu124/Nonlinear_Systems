function hoplotf(dir,pt,wnr,cnr,varargin)
% HOPLOTF: call hoplot for point read from dir/pt
%
%  hoplotf(dir,pt,wnr,cnr,varargin)
%
% varargin passed through to hopot
try p=loadp(dir,pt); 
catch fprintf(['point ' dir '/' pt ' does not exist.\n']); 
    npt=ptselect(dir); p=loadp(dir,npt); pt=npt; 
end
if p.sw.para==6; T=p.u(p.nu+p.hopf.iT); else T=p.hopf.T; end 
fprintf('lam=%g, T=%g\n',getlam(p),T); 
if nargin>4; hoplot(p,wnr,cnr,varargin{1});  
else try; hoplot(p,wnr,cnr,p.hopf.aux); catch; hoplot(p,wnr,cnr);  end 
end
