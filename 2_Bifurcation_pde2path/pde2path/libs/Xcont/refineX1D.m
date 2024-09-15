function [p,flag]=refineX1D(p,varargin)
% refineX1D: mesh adaptation, for 1D mfds (curves) 
%  p=refineX(p,varargin)
ngen=p.nc.ngen; imaxs=p.nc.imax; p.nc.imax=2*p.nc.imax; 
try maxt=p.nc.maxt; catch; maxt=2*size(p.tri,1); end 
sig=p.nc.sig;  if nargin>1; sig=varargin{1}; end
onp=p.np; flag=0;  
for j=1:ngen  % loop over refinements 
  u=p.u; par=u(p.nu+1:end); lamd=p.tau(p.nu+1:end);   
  idx=p.fuha.e2rs(p,sig); % select elements to refine 
  fprintf('refining %i elements\n',length(idx)); 
  [p,un,taun]=Xref1D(p,idx); p.sol.xi=1/p.np;  
  p.nu=p.nc.neq*p.np; p.u=[un;par]; 
  %p.X=X; p.tri=trin; %interp1(oldp(:),p.X,npo,'pchip');                
  if size(lamd,1)>size(lamd,2); lamd=lamd'; end 
  p.tau=[taun;lamd'];   %size(p.u), size(p.tau), size(lamd)
  pplot1D(p,11); 
  r=resi(p,p.u); fprintf('inires=%g\n',norm(r,Inf));   % was: r=p.r
  [u1,r,iter,Gu,Glam,p]=nloop(p,p.u); cres=norm(r,p.sw.norm); fprintf('res=%g\n',cres);  
  if cres<p.nc.tol; flag=1; pplot(p,12); 
      [p,u]=updX(p,u1); end % refinement successful 
  if p.nt>maxt; fprintf('nt=%i>ntmax=%i, stopping adaption\n',p.nt,maxt);   break; end 
end
p.nc.imax=imaxs; fprintf('refined from %i to %i\n',onp,p.np);
