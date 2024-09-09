function [p,flag]=oomeshada(p,varargin)
% oomeshada: special version for 1D qfem
ngen=p.nc.ngen; noa=nargin-2; k=1;  imaxs=p.nc.imax; p.nc.imax=2*p.nc.imax; 
try maxt=p.nc.maxt; catch maxt=2*size(p.pdeo.grid.t,2); end 
try gosafe=p.sw.gosafe; catch; gosafe=0; end % safety-switch: if gosafe=1, successful 
% refinements are saved to disk to later reload the last successful one
% (necessary cause p.pdeo is a handle object and hence cannot be easily copied)
try swapfn=p.file.swapfn; catch; swapfn='maswap'; end 
while k<=noa 
   switch lower(varargin{k})       
       case 'ngen'; ngen=varargin{k+1}; k=k+2;
       case 'sig';  p.nc.sig=varargin{k+1}; k=k+2;
       otherwise; break
   end
end
oonp=p.np; flag=0; 
if gosafe; save(swapfn,'p'); end % save to later recover successful steps 
for j=1:ngen  % loop over refinements 
  p=simplecoarselob(p);  
  u=p.u; uf=p.mat.fill*u(1:p.nu); oldp=p.pdeo.grid.p; onp=p.np; 
  par=u(p.nu+1:end); tau=p.tau(1:p.nu)'; tauf=p.mat.fill*tau; 
  lamd=p.tau(p.nu+1:end);     
  [p,idx]=p.fuha.e2rs(p,p.u); % select elements to refine 
  fprintf('refining %i elements\n',length(idx)); p.pdeo.grid.refineMesh(idx);   
   p=two2lob(p); npo=p.pdeo.grid.p; np=p.pdeo.grid.nPoints; p.sol.xi=1/np; 
  ui=zeros(p.nc.neq*np,1); taui=ui;  
  for i=1:p.nc.neq % interpolate u and tau to new mesh 
    switch size(oldp,1)  
      case 1; ui((i-1)*np+1:i*np)=interp1(oldp(:),uf((i-1)*onp+1:i*onp),npo);        
         taui((i-1)*np+1:i*np)=interp1(oldp(:),tauf((i-1)*onp+1:i*onp),p.pdeo.grid.p(:));          
    end
  end
  p.np=np; p.nu=p.nc.neq*p.np; p.u=[ui;par]; %size(ui), size(taui), size(lamd)
  if size(lamd,1)>size(lamd,2); lamd=lamd'; end 
  p=box2per(p); if any(isnan(p.mat.M)) isnan(p.mat.M); end 
  p.tau=[p.mat.drop*taui;lamd']; 
  if size(oldp,1)==1 && p.sw.verb>1 % display new point distribution (1D) 
    figure(6); clf; plot(npo,'r*'); axis tight; 
    figure(7); clf; plot(npo,p.u(1:p.np),'*'); 
  else figure(6); clf; plotsol(p,6,1,1); axis(p.plot.axis); 
  end
  r=resi(p,p.u); 
  fprintf('inires=%g\n',norm(r,Inf));   % was: r=p.r
  [p.u,r]=nloop(p,p.u); cres=norm(r,p.sw.norm); fprintf('res=%g\n',cres);
  if cres<p.nc.tol; flag=1; if gosafe; save(swapfn,'p'); end; end % refinement successful 
  nt=size(p.pdeo.grid.t,2); 
  if nt>maxt; fprintf('nt=%i>ntmax=%i, stopping adaption\n',nt,maxt);   break; end 
end
p.xo=p.pdeo.grid.p; 
if gosafe; s=load(swapfn,'p'); p=s.p; end % recover last successful step 
p.nc.imax=imaxs; fprintf('refined from %i to %i\n',oonp,p.np);
