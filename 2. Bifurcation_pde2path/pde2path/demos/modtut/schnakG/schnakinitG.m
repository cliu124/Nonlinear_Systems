function p=schnakinitG(p,par,np,sw) % init with prep. of mesh-adaption and usrlam
p=stanparam(p); screenlayout(p); p.nc.neq=2; p.sw.sfem=-1; 
G=mygraph(np,sw); p.G=G; p.np=G.numnodes; p.nu=2*p.np; p.sol.xi=1/(p.nu);
p=oosetfemops(p); lam=par(1); u=lam*ones(p.np,1); v=(1/lam)*ones(p.np,1); 
p.u=[u;v;par']; % initial solution guess with parameters
p.usrlam=0.5:0.5:4; % compute point and write to disk at these par values
p.nc.nsteps=20; p.sw.foldcheck=1; p.mesh.maxt=100; p.nc.neig=min(10,p.np); 
p.plot.auxdict={'\lambda','\sigma','d'}; 
p.nc.ilam=2; p.nc.lammin=-0.2; p.nc.lammax=1; p.plot.pstyle=-1; 
p.sol.ds=0.1; p.nc.dsmax=0.1; 