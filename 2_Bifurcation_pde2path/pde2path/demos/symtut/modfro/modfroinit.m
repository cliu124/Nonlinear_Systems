% initialization for modulated fronts after Balmforth et al 1999
function p=modfroinit(p,lx,nx,par,dir) 
p=stanparam(p); if ~exist(['./' dir],'dir'); mkdir(['./' dir]); end;  
p=setfn(p,dir); screenlayout(p);
p.nc.neq=2; p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
p.fuha.outfu=@hobra; vol=2*lx; 
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; p.vol=2*lx; 
p.np=pde.grid.nPoints; p.nu=2*p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e; p.mesh.nt=size(t,2);
p=setfemops(p); p.nc.ilam=1; p.nc.lammin=0; p.nc.lammax=1; p.nc.nsteps=20; 
x=getpte(p); x=x'; u=0.5*(1-tanh(2*x)); v=1-u; p.u=[u;v; par']; % initial guess 
p.u0=[u;v]; p.u0x=p.Dx*p.u0; % use same u as reference wave 
p.nc.sf=1; % stiffness for DBC (u,v)_=(1,0), (u,v)_+=(0,1) 
p.sol.ds=0.01; p.nc.dsmin=0.0001; p.nc.dsmax=0.05; p.sw.verb=2; p.sw.bifcheck=2; 
p.plot.pcmp=[1 2]; p.plot.pbcmp=3; p.plot.cl={'black','blue'};  plotsol(p); 