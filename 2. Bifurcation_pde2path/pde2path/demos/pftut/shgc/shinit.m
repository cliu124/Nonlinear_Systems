function p=shinit(p,nx,lx,ly,ndim,par,varargin)  % GCSH as 2 component system 
p=stanparam(p); p.nc.neq=2; p.ndim=ndim; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
p.fuha.lss=@gclss; p.fuha.blss=@gcblss; % Sherman-Morrison versions 
p.sw.runpar=0;  % switch off parfor in pmnewtonloop (clashes with global vars) 
p.sw.eigssol=1; % use Sherman-Morrison in eigs (gclsseigs) 
p.fuha.outfu=@shbra; p.sw.spcalc=0; p.plot.bpcmp=0; 
p.sol.ds=0.01; p.sol.dsmax=0.1; p.pm.resfac=1e-3; p.sw.bifcheck=2; 
p.sw.spcont=0; p.sw.spcalc=1; p.sw.sfem=-1; 
switch ndim 
  case 1; pde=stanpdeo1D(lx,2*lx/nx); p.Om=2*lx; p.fuha.outfu=@shbra;
  case 2; ny=round(nx*ly/lx); try sw=varargin{1}; catch; sw=[]; end; 
      pde=stanpdeo2D(lx,ly,nx,ny,sw); p.Om=4*lx*ly; p.plot.pstyle=2; 
      p.plot.axis='image'; p.sol.ds=0.001; p.sol.dsmax=0.01; p.sol.dsmin=0.001; 
  case 3; lz=varargin{1}; try sw=varargin{2}; catch; sw=[]; end 
      pde=stanpdeo3D(lx,ly,lz,nx,round(nx*ly/lx),round(nx*lz/lx),sw);
      p.Om=8*lx*ly*lz; p.plot.pstyle=1; p.plot.axis='image'; 
end
p.np=pde.grid.nPoints;  p.pdeo=pde; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;
p.nc.lammin=-4; p.nc.lammax=2; p.nc.ds=0.01; p.nc.dsmax=0.1; 
u=0*ones(p.np,1); v=u; u0=[u v]; p.u=u0(:); 
p.u=[p.u; par]; p.nc.ilam=1; p.file.smod=5; p.nc.neig=20;
screenlayout(p); p=setfemops(p); p.nc.nsteps=100;