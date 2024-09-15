function p=chinit(p,lx,nx,par) % Cahn-Hilliard
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; nx=nx(1); 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.sw.jac=1; p.fuha.outfu=@Jbra; 
switch size(lx,2)
    case 1; pde=stanpdeo1D(lx,2*lx/nx); p.dim=1; p.Om=2*lx; 
    case 2; sw.sym=2; pde=stanpdeo2D(lx(1),lx(2),nx,nx,sw); p.plot.pstyle=2; p.dim=2; 
      p.Om=4*lx(1)*lx(2); p.plot.axis='image'; 
  case 3; sw.sym=0; pde=stanpdeo3D(lx(1),lx(2),lx(3),0.1); %nx,nx,nx,sw); 
    p.plot.pstyle=2; p.dim=3; 
      p.Om=8*lx(1)*lx(2)*lx(3); p.plot.axis='image';  p.plot.pstyle=3; 
end
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; 
p.sol.xi=1/(p.nu); p.u=par(1)*ones(p.np,1); p.u=[p.u; par']; % initial sol, with param appended 
p=oosetfemops(p); p.nc.nsteps=20; p.sw.foldcheck=1; p.plot.axis='image'; 
p.plot.auxdict={'m','eps','lam'}; % parameter names 
p.plot.bpcmp=5; p.nc.lammin=-1.2; p.nc.lammax=1.2; p.nc.mu1=1; p.sw.bifcheck=2; 
p.sol.ds=0.1; p.nc.dsmax=0.05; p.sw.bprint=3; p.nc.nsteps=100; p.nc.neig=30;
p.nc.nq=1; p.fuha.qf=@qf; p.fuha.qfder=@qfder; p.sw.qjac=1; p.nc.ilam=[1,3]; % aux eqns 
p.fuha.spjac=@spjac;  p.sw.spjac=1; % for spectral cont 
p0.usrlam=[-0.5 0 0.5]; 