function p=FHNinit(p,lx,nx) % initialization routine for fhneps demo
% initialize structure, directory and screen:
p=stanparam(p); dir='init'; if ~exist(['./' dir],'dir'); mkdir(['./' dir]); end;  
p=setfn(p,dir); screenlayout(p);
% PDE settings:
p.nc.neq=2; p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
% boundary and mesh:
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; 
p.np=pde.grid.nPoints; p.nu=2*p.np; p.sol.xi=1/(p.nu); 
% background mesh for refinement:
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e; p.mesh.nt=size(t,2);
p=setfemops(p); 
% initial settings for continuation:
p.nc.ilam=1; p.nc.lammin=0; p.nc.lammax=10; p.nc.nsteps=20; 
p.sol.ds=0.01; p.nc.dsmin=0.0001; p.nc.dsmax=0.1; p.nc.imax=40; 
% plotting:
p.plot.auxdict={'eps','vel','p_3','p_4','p_5','p_6'}; 
p.plot.pcmp=[1 2]; p.plot.cl={'black','blue'};