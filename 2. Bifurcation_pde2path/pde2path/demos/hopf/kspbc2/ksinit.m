function p=ksinit(p,nx,lx,par) 
p=stanparam(p); p.nc.neq=2; p.ndim=1; 
p.fuha.sG=@sG;p.fuha.sGjac=@sGjac; 
p.fuha.outfu=@hobra; p.sw.spcalc=0; p.plot.bpcmp=0; 
p.nc.nq=1; p.fuha.qf=@qf; p.fuha.qfder=@qjac;
pde=stanpdeo1D(lx,2*lx/nx); p.np=pde.grid.nPoints; p.pdeo=pde;
p.vol=2*lx; p.sw.sfem=-1; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;
p.nc.lammin=1e-4;  p.sw.spcalc=1;  p.sw.bifcheck=2; p.nc.bisecmax=20; 
u=par(2)*ones(p.np,1); v=u; u0=[u v]; p.u=u0(:); 
p.u=[p.u; par]; p.file.smod=0; p.sol.ds=-0.01; p.nc.dsmax=0.02;
p=box2per(p,1); 



