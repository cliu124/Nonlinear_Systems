function p=schnakinit(p,dom,nx,par,varargin)
p=stanparam(p); screenlayout(p); p.nc.neq=2; p.sw.sfem=-1;
p.plot.auxdict={'\lambda','c','d','m'}; p.fuha.outfu=@hobra;
switch length(dom)
  case 1; lx=dom; p.pdeo=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; 
  case 2; ny=round(dom(2)/dom(1)*nx); lx=dom(1); ly=dom(2); p.vol=4*lx*ly; 
      pde=stanpdeo2D(lx,ly,nx,ny,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
end
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; p.nc.neig=10; p.nc.nsteps=50; 
p=setfemops(p); p.sol.xi=1/p.nu; p.file.smod=5; p.sw.para=2; p.sw.foldcheck=1; 
p.nc.ilam=1; p.sol.ds=-0.1; p.nc.dsmax=0.5; p.nc.lammin=0.1; p.sw.bifcheck=2; 
lam=par(1); m=par(4); al=m; bet=1-2*m; 
u=lam^al*ones(p.np,1); v=lam^bet*ones(p.np,1); 
p.u=[u;v;par']; % initial solution guess with parameters
[po,tr,ed]=getpte(p); p.mesh.bp=po; p.mesh.be=ed; p.mesh.bt=tr; p.nc.ngen=1; 
p.plot.bpcmp=6; p.plot.axis='image'; p.plot.cm='parula'; 
p.nc.resfac=1e-3; p.pm.mst=4; p.nc.del=1e-8; % small delta helpful for FDs 