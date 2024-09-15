function p=acinit(p,R,phi,nr,par) % AC on sector 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.sw.jac=1;
pde=secpdeo(R,phi,nr); p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; 
p.plot.pstyle=2; p.plot.cm='cool'; p.plot.axis='image'; p.nc.neig=20; 
p.sol.xi=1/(p.nu); p.u=zeros(p.np,1); p.u=[p.u; par']; % initial sol, with param appended 
plotsol(p); p=oosetfemops(p); p.nc.nsteps=20; p.sw.foldcheck=1; 
p.plot.auxdict={'lambda','s','gamma'}; % parameter names 
p.nc.lammax=4; p.sol.ds=0.05; p.nc.dsmax=0.2; p.nc.nsteps=100; 
p.usrlam=[0 2 4]; p.vol=2*phi*R^2; 