function p=schnakinit(p,nper,ppp,par)
p=stanparam(p); screenlayout(p); p=setfn(p,'p'); p.nc.lammin=2.7; 
p.nc.neq=2; p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
kc=sqrt(sqrt(2)-1); lx=nper*pi/kc; % wavenumer of the critical mode
p.pdeo=stanpdeo1D(lx,2*pi/kc/ppp); p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; 
p=setfemops(p); p.nc.ilam=1; p.sol.xi=1/p.nu; p.sol.ds=-0.01; p.nc.dsmax=0.1; 
u=par(1)*ones(p.np,1); v=(1/par(1))*ones(p.np,1); % hom.soln 
p.u=[u;v;par']; p.nc.nsteps=100; p.plot.pmod=0; p.file.smod=0; 