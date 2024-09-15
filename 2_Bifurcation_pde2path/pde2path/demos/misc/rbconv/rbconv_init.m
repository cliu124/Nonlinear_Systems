function p=rbconv_init(p) 
% initialization routine for "rbconv" demo
% see manual or routine "stanparam" for explanations of settings

% standard settings, window arrangement, file names, plotting
p=stanparam(p); dir=sprintf('%s',inputname(1)); p=setfn(p,dir); screenlayout(p);
p.plot.pstyle=2; p.plot.pcomp=2; 

% PDE standard implementation
p.nc.neq=3; p.fuha.G=@rbconv_G; p.fuha.Gjac=@rbconv_Gjac; 
p.fuha.bc=@rbconv_bc_noslip; p.fuha.bcjac=@rbconv_bc_noslip;

% domain and mesh
lx=2; lz=0.5; p.mesh.geo=rec(lx,lz); 
p=stanmesh(p,100,25); p=setbmesh(p); 

% continuation settings
p.nc.ilam=1;
p.sol.xi=1/p.np; p.sol.ds=0.1; p.nc.dsmin=0.00001; p.nc.dlammax=50; p.nc.nsteps=10; p.nc.neig=20; 
p.nc.amod=0; p.nc.tol=1e-10; 
p.nc.lammax=900; 

% initial point and parameters
p.u=0*ones(p.nc.neq*p.np,1); 
par(1)=740; 
p.u = [p.u; par'];

p=setfemops(p);