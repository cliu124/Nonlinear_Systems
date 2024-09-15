function p=gpinit(p,lx,nx,par) % init 1D GP, real scalar version 
p=stanparam(p); screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
p.fuha.outfu=@gpbra;
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; % domain and mesh
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;
p.u=zeros(p.np,1); p.u=[p.u; par']; 
p=setfemops(p); p.nc.nsteps=20; p.sw.foldcheck=1; 
p.plot.auxdict={'lam','s','sig','ga'}; 
p.plot.pstyle=1; p.plot.cm='cool'; p.usrlam=-1:1:5; 
p.nc.nsteps=100; p.sw.bifcheck=2; p.sw.jac=1; 
p.nc.ilam=1; p.nc.lammax=5; p.sol.ds=0.1; p.nc.dsmax=0.2;
p.nc.maxt=600; p.nc.sig=0.6; p.nc.ngen=4; % meshada settings 
p.usrlam=0:1:5; p.nc.mu1=0.5; p.nc.sf=0; 