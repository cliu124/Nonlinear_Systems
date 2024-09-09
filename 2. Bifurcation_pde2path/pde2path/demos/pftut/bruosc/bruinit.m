function p=bruinit(p,dom,nx,par,varargin)
p=stanparam(p); screenlayout(p); p.nc.neq=2; p.sw.sfem=-1; 
p.plot.auxdict={'a','b','D_u','D_v'}; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.spjac=@spjac; p.fuha.outfu=@hobra; 
switch length(dom)
  case 1; lx=dom; p.pdeo=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; p.file.smod=5; 
  case 2; lx=dom(1); ly=dom(2); ny=round(ly/lx*nx); p.vol=4*lx*ly; p.file.smod=5; 
      pde=stanpdeo2D(lx,ly,nx,ny,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
end
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; p.nc.neig=30; p.nc.nsteps=50; 
p=setfemops(p); p.sol.xi=1/p.nu; p.sw.para=2; p.sw.foldcheck=1; 
p.nc.ilam=3; p.sol.ds=-0.1; p.nc.dsmax=0.5; p.nc.lammin=2; p.sw.bifcheck=2; 
a=par(1); b=par(2); u=a*ones(p.np,1); v=(b/a)*ones(p.np,1); 
p.u=[u;v;par']; % initial solution guess with parameters
[po,tr,ed]=getpte(p); p.mesh.bp=po; p.mesh.be=ed; p.mesh.bt=tr; p.nc.ngen=1; 
p.plot.bpcmp=0; p.plot.axis='image'; p.plot.cm='hot'; p.plot.fancybd=2; 
p.nc.resfac=1e-3; p.pm.mst=4; 