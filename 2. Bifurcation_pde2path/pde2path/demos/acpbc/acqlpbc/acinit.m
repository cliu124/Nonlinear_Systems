function p=acinit(p,lx,nx,par,dim,per) % acql 
p=stanparam(p); p.nc.neq=1; p.nc.ilam=2; screenlayout(p); 
p.sw.sfem=-1; p.dim=dim; p.fuha.sG=@sG; p.sw.jac=1; p.sw.bcper=per;
switch dim % domain setup 
 case 1; pde=stanpdeo1D(lx,2/nx); p.fuha.sGjac=@sGjac1D; 
 case 2; ly=0.9*lx; ny=round(nx*ly/lx); p.plot.pstyle=3; p.fuha.sGjac=@sGjac2D; 
     sw.sym=2; pde=stanpdeo2D(lx,ly,nx,ny,sw); % init with symmetric mesh
 case 3; ly=0.75*lx;lz=0.5*lx; pde=stanpdeo3D(lx,ly,lz,2*lx/nx); 
     p.fuha.sGjac=@sGjac3D; 
end
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); p.nt=pde.grid.nElements; 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;  % background mesh  
p.fuha.e2rs=@e2rs; % and refinement selector 
p.u=zeros(p.np,1); p.u=[p.u; par]; p.sw.foldcheck=0; p.sw.bifcheck=2; 
if per~=0; p=box2per(p,per); end
p.usrlam=1:1:4; p.nc.lammax=5;  p.sol.ds=0.2; p.nc.dsmax=0.4; p.nc.nsteps=20; 