function p=shinitpsq(p,nbd,del,hmax,par,varargin) % SH on perturbed square
p=stanparam(p); p.nc.neq=2; p.nc.neig=20; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
p.sw.spcalc=0; p.plot.bpcmp=0; 
p.sol.ds=0.01; p.sol.dsmax=0.1; p.pm.resfac=1e-3; p.sw.bifcheck=2; 
p.sw.spcont=0; p.sw.spcalc=1; p.sw.sfem=-1; p.sw.foldcheck=0; 
t=linspace(0,2*pi,nbd); t=t(1:end-1); 
x1=t; y1=del*sin(t/2); x3=2*pi-t; y3=2*pi-del*sin(t/2); 
x2=2*pi-del*sin(t/2); y2=t;  x4=del*sin(t/2); y4=2*pi-t; 
x=[x1,x2,x3,x4]; y=[y1,y2,y3,y4];
pde=freegeompdeo(x,y,hmax); % h as argument
p.plot.axis='image'; p.sol.ds=0.001; p.sol.dsmax=0.01; p.sol.dsmin=0.001; 
p.np=pde.grid.nPoints;  p.pdeo=pde; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;
p.nc.lammin=-4; p.nc.lammax=2; p.nc.ds=0.01; p.nc.dsmax=0.1; 
u=0*ones(p.np,1); v=u; u0=[u v]; p.u=u0(:); 
p.u=[p.u; par]; p.nc.ilam=1; p.file.smod=5; p.plot.pstyle=2; 
%screenlayout(p); 
p=setfemops(p); p.nc.nsteps=100;