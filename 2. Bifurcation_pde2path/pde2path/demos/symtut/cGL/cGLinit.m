function p=cGLinit(p,lx,nx) % (generic) init routine for cGL problem 
p=stanparam(p); screenlayout(p); % set standard parameters and screenlayout 
p.nc.neq=2; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; % OOPDE setting of hom. Neumann BC 
p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;  
% background-mesh (for oomeshadac):
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e; p.mesh.nt=size(t,2); 
p=setfemops(p); p.nc.ilam=1; p.sw.foldcheck=1; p.sw.bifcheck=2; 
p.nc.neig=10; % very low due to poor convergence of evals
p.nc.mu1=0.5; % be relaxed about possible bif-detection
p.nc.nsteps=20; p.usrlam=0:0.5:2; p.nc.lammax=2;  % user-vals for output
p.plot.auxdict={'r','vel','\nu','\mu','c3','c5','\gamma','sym','scale'}; 
p.plot.pcmp=[1 2]; p.plot.cl={'black','blue'};
p.file.smod=10; p.sol.ds=0.1; p.nc.dsmax=0.5; % saving, stepsize, max stepsize 