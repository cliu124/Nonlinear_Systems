function p=hexinit(p,hmax,par) 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.nc.neig=20; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
pde=hexpdeo(hmax); p.pdeo=pde; % domain and mesh
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;
p.u=zeros(p.np,1); p.u=[p.u; par']; p=setfemops(p);   
p.nc.nsteps=20; p.sw.foldcheck=0; p.nc.ilam=1; p.nc.mu1=2; 
p.nc.lammax=50; p.sol.ds=0.1; p.nc.dsmax=0.5; p.sw.bifcheck=2;
p.sw.verb=2; p.plot.auxdict={'\lambda'}; 
p.plot.pstyle=2; p.plot.cm='cool'; p.plot.axis='image'; 