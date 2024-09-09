function p=acinit(p,lx,ly,nx,par) 
p=stanparam(p); screenlayout(p); p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
pde=stanpdeo2D(lx,ly,2*lx/nx); % % domain and mesh, h as argument
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
p.u=zeros(p.np,1); p.u=[p.u; par']; p.nc.nsteps=20; 
p=box2per(p,2); % prepare fill, drop for periodic BC, here in y
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;
p.plot.auxdict={'c','lambda','gamma','d'}; p.sw.foldcheck=1; 
p.plot.pstyle=1; p.plot.cm='cool'; p.usrlam=[0 0.5 1]; 
p.nc.nsteps=100; p.sw.bifcheck=2; p.sw.jac=1; 