function p=schnakcinit(p,nx,par,ell) % init on cone 
p=stanparam(p); screenlayout(p); p.nc.neq=2; p.sw.sfem=-1;
p.plot.auxdict={'\lambda','\sigma','D','a','eps', 's'}; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
pde=diskpdeo2(1,nx,4,ell); % ellipse-PDE object 
p.pdeo=pde; p.plot.pstyle=2; p.sw.verb=2; 
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; p.nc.neig=40; p.nc.nsteps=50; 
p.sol.xi=1/p.nu; p.file.smod=5; p.sw.para=2; p.sw.foldcheck=1; p.nc.dsmin=0.01; 
p.nc.ilam=1; p.sol.ds=-0.1; p.nc.dsmax=0.1; p.nc.lammin=2; p.sw.bifcheck=2; 
lam=par(1); u=lam*ones(p.np,1); v=(1/lam)*ones(p.np,1); p.u=[u;v;par']; 
p=oosetfemops(p); p.tau(1:p.nu+1)=[0*[u; u]; 0]; 
p.plot.bpcmp=7; p.plot.axis='image'; p.plot.cm='cool'; p.usrlam=2:0.1:3.3; 
p.pm.resfac=1e-4; p.nc.intol=-1e-5; 