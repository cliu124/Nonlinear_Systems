function p=bruinit(p,dom,nx,par,varargin)
p=stanparam(p); screenlayout(p); p.nc.neq=2; p.sw.sfem=-1; 
p.plot.auxdict={'a','b','D_u','D_v','\delta','\omega'}; 
p.fuha.sG=@sGpf; p.fuha.sGjac=@sGjacpf; p.fuha.outfu=@hobra; 
switch length(dom)
  case 1; lx=dom; p.pdeo=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; 
  case 2; lx=dom(1); ly=dom(2); ny=round(ly/lx*nx); p.vol=4*lx*ly; 
      pde=stanpdeo2D(lx,ly,nx,ny,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
end
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; p.nc.neig=30; p.nc.nsteps=50; 
p=setfemops(p); p.sol.xi=1/p.nu; p.sw.para=2; p.sw.foldcheck=1; p.file.smod=5; 
p.nc.ilam=3; p.sol.ds=-0.1; p.nc.dsmax=0.5; p.nc.lammin=-1; p.sw.bifcheck=2; 
a=par(1); b=par(2); u=a*ones(p.np,1); v=(b/a)*ones(p.np,1); 
p.nu=p.nu+2; p.u=[u;v;0;0;par']; % initial solution augmented by oscillator vars v
p.plot.bpcmp=0; p.plot.axis='image'; p.plot.cm='hot'; p.plot.shading='interp'; 