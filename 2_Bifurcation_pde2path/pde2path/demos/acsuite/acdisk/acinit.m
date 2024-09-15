function p=acinit(p,r1,nr,par) % init AC with OOPPDE 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.sw.jac=1;
%pde=diskpdeo(r1,nr);
pde=diskpdeo2(r1,20,4); % alternate simple version, all bd-segm=1
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; 
p.plot.pstyle=2; p.plot.cm='cool'; p.nc.neig=40; 
p.sol.xi=1/(p.nu); p.u=zeros(p.np,1); p.u=[p.u; par']; % initial sol, with param appended 
plotsol(p); p=oosetfemops(p); p.nc.nsteps=20; p.sw.foldcheck=1; 
p.plot.auxdict={'lambda','s','gamma'}; % parameter names 
p.nc.lammax=4; p.sol.ds=0.05; p.nc.dsmax=0.05; p.nc.nsteps=100; 
p.usrlam=[-1 0 1 2]; 