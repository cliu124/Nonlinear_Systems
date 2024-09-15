function p=acinit(p,lx,ly,nx,par) % ac on torus
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; p.sw.jac=1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
p.fuha.qf=@qf; p.fuha.qfder=@qfder; p.nc.xiq=0.1; 
sw.sym=2; pde=stanpdeo2D(lx,ly,nx,round(nx*ly/lx),sw); 
p.plot.auxdict={'R','rho','c','lambda','gamma','s'}; p.plot.pstyle=2; p.plot.cm='cool'; 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e; % for mesh-ref 
p.u=zeros(p.np,1); p.u=[p.u; par']; p.nc.ilam=3; p.nc.nsteps=20; %p.sw.foldcheck=1; 
p=box2per(p,[1 2]);  % switch on pBC in x and y 
p=setfemops(p); p.sol.ds=0.1; p.nc.dsmax=0.2; p.nc.lammax=2.2; p.usrlam=[2 5]; 