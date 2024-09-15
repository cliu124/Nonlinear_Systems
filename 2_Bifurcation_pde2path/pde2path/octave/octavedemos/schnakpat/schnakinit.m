function p=schnakinit(p,dom,nx,par,varargin)
p=stanparam(p); screenlayout(p); p.nc.neq=2; p.sw.sfem=-1;
p.sw.spjac=1; % use analytical Jacobian for spectral point cont (fold cont)
p.plot.auxdict={'\lambda','\sigma','d','||u_1||_{\infty}','min(|u_1|)'};
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.spjac=@spjac; 
try sw=varargin{1}; sym=sw.sym; catch; sym=0; end; sym 
switch length(dom)
  case 1; lx=dom; p.pdeo=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; 
  case 2; ny=round(dom(2)/dom(1)*nx); lx=dom(1); ly=dom(2); p.vol=4*lx*ly; 
      pde=stanpdeo2D(lx,ly,lx/nx,sym); p.pdeo=pde; p.plot.pstyle=2; 
  case 3; lx=dom(1); ly=dom(2); lz=dom(3); h=2*lx/(nx-1);  p.vol=8*lx*ly*lz; h
      pde=stanpdeo3D(lx,ly,lz,h,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
      p.plot.EdgeColor='none'; 
end
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; p.nc.neig=30; p.nc.nsteps=50; 
p=setfemops(p); p.sol.xi=1/p.nu; p.file.smod=10; p.sw.para=2; p.sw.foldcheck=1; 
p.nc.ilam=1; p.sol.ds=-0.1; p.nc.dsmax=0.1; p.nc.lammin=2; p.sw.bifcheck=2; 
lam=par(1); u=lam*ones(p.np,1); v=(1/lam)*ones(p.np,1); 
p.u=[u;v;par']; % initial solution guess with parameters
[po,tr,ed]=getpte(p); p.mesh.bp=po; p.mesh.be=ed; p.mesh.bt=tr; p.nc.ngen=1; 
p.plot.bpcmp=4; p.plot.axis='image'; p.plot.cm='hot'; p.plot.fancybd=2; 
p.nc.resfac=1e-3; p.pm.mst=4; 