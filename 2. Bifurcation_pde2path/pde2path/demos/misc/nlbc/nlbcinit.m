function p=nlbcinit(p,nx) % initialization routine for "nlbc" demo
% see manual or routine "stanparam" for explanations of settings
% standard settings, window arrangement, file names
p=stanparam(p); dir=sprintf('%s',inputname(1)); p=setfn(p,dir); screenlayout(p); 
% PDE standard implementation
p.nc.neq=1;          % set number of equations for pde
p.fuha.G=@G;         % set functions G in pde
p.fuha.Gjac=@Gjac;   % set jacobian of G (not required if p.jsw=3)
p.fuha.bc=@nlbc; p.fuha.bcjac=@nlbcjac; 
% PDE simple implementation
p.sw.sfem=0; p.eqn.c=1; p.eqn.b=0; p.eqn.a=0;
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
%p.mesh.geo=rec(1,1);ny=nx; p=stanmesh(p,nx,ny); % square interesting too
p.mesh.geo=circgeo(1,nx); hmax=2/nx; p=stanmesh(p,hmax);  % domain and mesh
p=setbmesh(p);
% continuation settings
p.nc.ilam=1;       % set primary parameter index (here no other active parameters)
p.sol.ds=0.5; p.nc.nsteps=10; p.nc.lammax=10;
p.nc.amod=0; p.sol.xi=1/p.np; 
% initial point and parameters
p.eqn.a=0.0; 
p.u=0*ones(p.nu,1); 
par(1)=0.1; p.u = [p.u; par'];