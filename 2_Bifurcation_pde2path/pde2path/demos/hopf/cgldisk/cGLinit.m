function p=cGLinit(p,lx,nx) % (generic) init routine for cGL problem 
p=stanparam(p); screenlayout(p); % set standard parameters and screenlayout 
p.nc.neq=2; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; p.fuha.outfu=@hobra; 
pde=diskpdeo(lx,nx); p.vol=2*pi*lx^2; p.ndim=2;
p.plot.pstyle=3; p.plot.viev=[30,50]; 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;  
p=setfemops(p); p.nc.ilam=1; p.sw.foldcheck=0; p.sw.bifcheck=2; 
p.nc.neig=20; p.nc.mu1=0.5; % be relaxed about possible bif-detection
p.nc.nsteps=20; p.usrlam=0:2:6; p.nc.lammax=8;  % user-vals for output
p.plot.auxdict={'r','\nu','\mu','c3','c5','s','\delta'}; 
p.plot.bpcmp=11; p.file.smod=10; p.sol.ds=0.1; p.nc.dsmax=0.5; % saving, stepsize, max stepsize 