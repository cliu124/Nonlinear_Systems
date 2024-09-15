function p=acinit(p,lx,ly,nx,ny,par,ref) % AC on sphere
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
p.fuha.qf=@qf; p.fuha.qfder=@qfder; 
pde=spherepdeo(lx,ly,nx,ny,ref,2); % spherical coord. pdeo; uses modified mesh 
p.plot.auxdict={'c','lambda','gamma','s'}; p.plot.pstyle=2; p.plot.cm='cool'; 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;
p.u=zeros(p.np,1); p.u=[p.u; par']; p.nc.ilam=2; p.nc.neig=40; p.usrlam=0:1:3;  
p=box2per(p,1); % switch on pBC in x
p.nc.nsteps=20; p.sw.foldcheck=1; p.sw.bifcheck=2; p.nc.mu1=0.5; p.nc.mu2=0.01; 
p.nc.ilam=2; p.nc.lammax=3.1; p.sol.ds=0.1; p.nc.dsmax=0.21; p.sw.verb=2; 