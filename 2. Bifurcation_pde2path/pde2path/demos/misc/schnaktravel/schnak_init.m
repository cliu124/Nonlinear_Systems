function p=schnak_init(p,mx,nx) 
% initialization routine for "schnaktravel" demo, OOPDE 
p=stanparam(p); dir=sprintf('%s',inputname(1)); p=setfn(p,dir); screenlayout(p);
p.nc.neq=2; p.sw.sfem=-1; 
p.fuha.sG=@schnak_sG; p.fuha.sGjac=@schnak_sGjac; p.fuha.e2rs=@e2rs;  

% domain and mesh
kc=sqrt(sqrt(2)-1);lx=mx*2*pi/kc/2; p.Om=2*lx; 
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; % domain and mesh
p.np=pde.grid.nPoints; p.nu=2*p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t;  
p.mesh.be=e; p.mesh.nt=size(t,2); % background-mesh (for oomeshadac)
p.nc.ilam=1; p.nc.lammin=2; p.nc.lammax=4; % continuation settings

lam=sqrt(60)*sqrt(3-sqrt(8))+0.1; % initial guess and parameters
u=lam*ones(p.np,1); v=(1/lam)*ones(p.np,1); p.u=[u;v]; p.sol.xi=1/(p.nu); 
par(1)=lam; % production parameter
par(2)=0;   % comoving frame
par(3)=1;  % diff.coeff 
p.u = [p.u; par']; p=setfemops(p); 
