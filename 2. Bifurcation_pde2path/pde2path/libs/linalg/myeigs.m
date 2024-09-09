function [V,mu]=myeigs(Gu,M,neig,om,opts,p)
% myeigs: interface for eigs, controlled by p.sw.eigssol 
% neig=p.nc.neig; om=p.nc.eigref(1); 
n=size(M,1); 
try eigssol=p.sw.eigssol; catch; eigssol=0; end 
switch eigssol
case 0 % STANDARD 
 if p.sw.eigsstart==1; p.sw.evopts.v0=ones(n,1)/n;  
 else; try; p.sw.evopts=rmfield(p.sw.evopts,'v0'); catch; end; end 
 [V,mu]=eigs(Gu,M,neig,om,opts);    
case 1 % global coupling version with global vars 
  opts.tol=1e-4; opts.p=2*p.nc.neig+1; n=size(M,1); 
  [V,mu]=eigs(@(b) gcafun(p,Gu,M,om,b),n,M,neig,om,opts); 
case 2 % ilu 
  opts.tol=1e-4; % opts.p=2*p.nc.neig+1; n=size(M,1); 
  [V,mu]=eigs(@(b) gcafunilu(p,Gu,M,om,b),n,M,neig,om,opts); 
case 3 % user-provided 
  opts.tol=1e-4; % opts.p=2*p.nc.neig+1; n=size(M,1); 
  [V,mu]=eigs(@(b) myeigsfu(p,Gu,M,om,b),n,M,neig,om,opts); 
end    