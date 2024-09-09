function p=acinit(p,lx,ly,nx,par) % ac2D 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
sw.sym=2; %p.pde=stanpdeo2D(lx,ly,2*lx/nx); % h as argument
pde=stanpdeo2D(lx,ly,nx,nx,sw); % alternate syntax with nx and ny
p.cl=0; % (convenience) switch for later comparison with classical FEM 
p.plot.auxdict={'c','lambda'}; p.plot.pstyle=2; p.plot.cm='cool'; p.plot.axis='equal'; 
p.pdeo=pde; % domain and mesh
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;
p.u=zeros(p.np,1); p.u=[p.u; par']; p.nc.ilam=2; p=setfemops(p);   
p.nc.nsteps=20; p.sw.foldcheck=1; p.sol.ds=0.1; p.nc.dsmax=0.2; 
p.usrlam=0:1:3; 
