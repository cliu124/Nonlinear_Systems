% init for breathers in FHN following Ikeda et al, Meth.Appl.Anal.2000
function p=brinit(p,lx,nx,par,pw) 
p=stanparam(p); dir='init'; if ~exist(['./' dir],'dir'); mkdir(['./' dir]); end;  
p=setfn(p,dir); screenlayout(p); 
p.nc.neq=2; p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
p.fuha.outfu=@hobra; 
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; p.vol=2*lx; 
p.np=pde.grid.nPoints; p.nu=2*p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e; p.mesh.nt=size(t,2);
x=getpte(p); x=x'; al=par(1); bet=par(2); ga=par(3);
um=(al+bet)/2-sqrt((al+bet)^2/4-(al*bet+1/ga)); 
up=(al+bet)/2+sqrt((al+bet)^2/4-(al*bet+1/ga)); 
u=um*(abs(x)>pw)+up*(abs(x)<pw+0.1); v=u/ga; 
skew=0; % use skew~=0 to converge to a trav.pulse 
p.u=[u+skew*x.*(abs(x)<=pw); v; par']; 
p=setfemops(p); 
p.nc.ilam=1; p.nc.lammin=0; p.nc.nsteps=20; p.sw.verb=2; p.sw.bifcheck=2; 
p.sol.ds=0.01; p.nc.dsmin=0.0001; p.nc.dsmax=0.04; p.nc.imax=10; 
p.plot.auxdict={'alpha','beta','ga','del','D','s'}; 
p.plot.pcmp=[1 2]; p.plot.bpcmp=8; p.plot.cl={'black','blue'};