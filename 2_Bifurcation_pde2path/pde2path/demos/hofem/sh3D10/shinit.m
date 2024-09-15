function p=shinit(p,lx,ly,lz,nx,ny,nz,par,hofem) % Swift-Hohenberg as 2 component system, with 
% singular M=diag(M,0) mass-matrix to compute correct eigenvalues 
p=stanparam(p); p.nc.neq=2; p.ndim=2; p.nc.neig=4; p.hofem=hofem; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
p.fuha.outfu=@shbra; p.sw.spcalc=0; p.plot.bpcmp=0; 
p.sol.ds=0.01; p.sol.dsmax=0.1; p.pm.resfac=1e-3; p.sw.bifcheck=2; 
p.sw.spcont=0; p.sw.spcalc=1; p.sw.sfem=-1; p.sw.foldcheck=1;  
% create stanpdeo3D, convert to 10node triangulation, 
sw.sym=1; pde=stanpdeo3D(lx,ly,lz,nx,ny,nz,sw); p.pdeo=pde; 
p.pdeo.grid.nPoints, 
p=four2ten(p); p.np, pause 
p.plot.pstyle=3; p.plot.cm='cool'; p.plot.axis='image'; p.sw10=1; 
p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;
p.nc.lammin=-4; p.nc.lammax=2; p.nc.ds=0.01; p.nc.dsmax=0.1; 
u=0*ones(p.np,1); v=u; u0=[u v]; p.u=u0(:); 
p.u=[p.u; par']; p.nc.ilam=1; p.file.smod=5; 
p=setfemops(p); p.nc.nsteps=100; 
p.Om=sum(sum(p.mat.M(1:p.np,1:p.np),1)); p.plot.shsw=0; 
plotsol(p,1,1,3); p.nc.mu2=0.005; figure(10); clf; spy(p.mat.Ks)