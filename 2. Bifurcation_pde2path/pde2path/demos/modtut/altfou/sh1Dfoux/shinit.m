function p=shinit(p,lx,nx,par) % init with prep. of mesh-adaption and usrlam
p=stanparam(p); screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
p.np=nx; p.lx=lx; p.nu=p.np; p.sol.xi=1/(p.nu);
p.u=zeros(p.np,1); p.u=[p.u; par']; p.usrlam=-0.5:0.5:1; 
p=oosetfemops(p);    % generate diff matrices 
p.nc.nsteps=20; p.sw.foldcheck=1; p.mesh.maxt=100; 
p.plot.auxdict={'\lambda','c2','c3'}; 
p.nc.ilam=1; p.nc.lammin=-1; p.nc.lammax=2; p.sol.ds=0.1; p.nc.dsmax=0.1; 