function p=cGLinit(p,lx,nx,par) % (generic) init routine for cGL problem 
p=stanparam(p); screenlayout(p); % set standard parameters and screenlayout 
p.nc.neq=2; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; p.fuha.outfu=@hobra; 
switch length(lx)
    case 1; pde=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; p.ndim=1; 
    case 2; sw.sym=0; pde=stanpdeo2D(lx(1),lx(2),nx,nx,sw); p.vol=4*lx(1)*lx(2); p.ndim=2;
end
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;  
p=setfemops(p); p.nc.ilam=1; p.sw.foldcheck=0; p.sw.bifcheck=2; 
p.nc.neig=20; p.nc.mu1=0.5; % be relaxed about possible bif-detection
p.nc.nsteps=20; p.usrlam=0:2:6; p.nc.lammax=8;  % user-vals for output
p.plot.auxdict={'r','\nu','\mu','c3','c5','s','\delta'}; 
p.plot.bpcmp=11; p.file.smod=10; p.sol.ds=0.1; p.nc.dsmax=0.5; % saving, stepsize, max stepsize 
u=zeros(p.np,1); v=u; p.u=[u;v; par]; % initial guess (here trivial) and pars 