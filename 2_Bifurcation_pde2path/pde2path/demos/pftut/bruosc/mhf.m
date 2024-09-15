function mhf(dir,pt,wnr,cnr,varargin)
% modified hoplotf
try p=loadp(dir,pt); 
catch fprintf(['point ' dir '/' pt ' does not exist.\n']); 
    npt=ptselect(dir); p=loadp(dir,npt); pt=npt; 
end
fprintf('lam=%g\n',getlam(p)); 
if nargin>4; myhoplot(p,wnr,cnr,varargin{1});  
else myhoplot(p,wnr,cnr);  end 
figure(1); colormap hot; colorbar; figure(6); 
