function [p,flag]=coarsenX1D(p,varargin)
% refineX1D: mesh adaptation, for 1D mfds (curves) 
%  p=refineX(p,varargin)
sig=p.nc.sig;  if nargin>1; sig=varargin{1}; end
onp=p.np; flag=0; ngen=p.nc.ngen; 
for j=1:ngen  % loop over refinements 
  u=p.u; par=u(p.nu+1:end); lamd=p.tau(p.nu+1:end);   
  idx=p.fuha.e2cs(p,sig); % select elements to refine   
  [p,un,taun]=Xcoarse1D(p,idx); p.sol.xi=1/p.np;  
  p.nu=p.nc.neq*p.np; p.u=[un;par]; 
  if size(lamd,1)>size(lamd,2); lamd=lamd'; end 
  p.tau=[taun;lamd'];   %size(p.u), size(p.tau), size(lamd)
  pplot1D(p,11); 
  r=resi(p,p.u); fprintf('inires=%g\n',norm(r,Inf));   % was: r=p.r
  [u1,r,iter,Gu,Glam,p]=nloop(p,p.u); cres=norm(r,p.sw.norm); fprintf('res=%g\n',cres);  
end
fprintf('coarsened from %i to %i\n',onp,p.np);
