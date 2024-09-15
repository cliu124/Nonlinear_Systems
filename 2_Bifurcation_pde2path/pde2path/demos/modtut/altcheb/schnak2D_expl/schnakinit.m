function p=schnakinit(p,lx,ly,nx,ny,par,Lfn)
p=stanparam(p); screenlayout(p); p.nc.neq=2; p.sw.sfem=-1; p.nc.neq=2; 
p.plot.auxdict={'\lambda','\sigma','d','||u_1||_{\infty}','min(|u_1|)'};
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.Lfn=Lfn; % filename for L 
p.nx=nx; p.ny=ny; p.np=nx*ny; p.lx=lx;  p.ly=ly; 
p.nu=p.nc.neq*p.np; p.sol.xi=1/(p.nu);
p.u=zeros(p.np,1); p.u=[p.u; par']; 
p.nu=p.np*p.nc.neq; p.nc.neig=30; p.nc.nsteps=50; 
p=oosetfemops(p); 
p.sol.xi=1/p.nu; p.file.smod=10; p.sw.para=2; p.sw.foldcheck=1; 
p.nc.ilam=1; p.sol.ds=-0.1; p.nc.dsmax=0.1; p.nc.lammin=2; p.sw.bifcheck=2; 
lam=par(1); u=lam*ones(p.np,1); v=(1/lam)*ones(p.np,1); 
p.u=[u;v;par']; % initial solution guess with parameters
p.plot.pstyle=-1; p.plot.bpcmp=4; p.plot.axis='image'; p.plot.cm='hot'; 
p.nc.resfac=1e-3; p.pm.mst=4; 