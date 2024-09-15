function hoplotfm(dir,pt,wnr,cnr,varargin)
try p=loadp(dir,pt); 
catch fprintf(['point ' dir '/' pt ' does not exist.\n']); 
    npt=ptselect(dir); p=loadp(dir,npt); pt=npt; 
end
fprintf('lam=%g\n',getlam(p)); 
if nargin>4; hoplot(p,wnr,cnr,varargin{1});  
else try; hoplot(p,wnr,cnr,p.hopf.aux); catch; hoplot(p,wnr,cnr);  end 
end
figure(wnr); view([0,90]); shading interp; axis image; 
title([dir '/' pt]); colorbar; ylabel('t'); 
