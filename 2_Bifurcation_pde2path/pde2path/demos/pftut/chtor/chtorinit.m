function p=chtorinit(p,lx,ly,nx,par)
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1;
p.plot.auxdict={'m','eps','lam', 'R','rho','s'}; % parameter names 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
p.fuha.qf=@qf2; p.fuha.qfder=@qf2der; p.fuha.outfu=@chbra; 
sw.sym=2; pde=stanpdeo2D(lx,ly,nx,round(nx*ly/lx),sw); 
p.pdeo=pde; p.plot.pstyle=2; 
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; p.nc.neig=20; 
p.sol.xi=1/p.nu; p.file.smod=10; p.sw.para=2; p.sw.foldcheck=1; 
p.sol.ds=0.1; p.nc.dsmax=0.1; p.nc.lammin=-1.2;  p.nc.lammax=1.2; 
p.sw.bifcheck=2; p.nc.lammax=0.1; p.nc.nsteps=100; 
p.u=[ones(p.np,1); par']; p=box2per(p,[1 2]); 
p.vol=p.mat.vMs*ones(p.nu,1); 
p.u=ones(p.nu,1)/p.vol; p.u=[p.u; par']; % initial sol, with param appended 
[po,tr,ed]=getpte(p); p.mesh.bp=po; p.mesh.be=ed; p.mesh.bt=tr; p.nc.ngen=1; 
p.plot.bpcmp=7; p.plot.axis='image'; p.plot.cm='cool'; p.usrlam=0; 