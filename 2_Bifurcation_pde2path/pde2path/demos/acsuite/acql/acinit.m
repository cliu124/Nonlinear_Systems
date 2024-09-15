function p=acinit(p,lx,nx,par,dim) % acql 
p=stanparam(p); p.nc.neq=1; p.nc.ilam=2; screenlayout(p); 
p.sw.sfem=-1; p.dim=dim; p.fuha.sG=@sG; p.sw.jac=1; 
p.jacsw=1; % additional switch for Jac: 1=approx., 
% 2=via numjac for \pa_u(div(c(u)\nab v), see sGjac*D for more comments
p.fuha.cfu=@cfu; % function handle for the nonlinear diffusion coeff
switch dim % domain setup 
 case 1; pde=stanpdeo1D(lx,2/nx); p.fuha.sGjac=@sGjac1D; 
 case 2; ly=0.9; ny=round(nx*ly/lx); p.plot.pstyle=3; p.fuha.sGjac=@sGjac2D; 
     sw.sym=2; pde=stanpdeo2D(lx,ly,nx,ny,sw); % init with symmetric mesh
 case 3; ly=0.75*lx;lz=0.5*lx; pde=stanpdeo3D(lx,ly,lz,2*lx/nx); 
     p.fuha.sGjac=@sGjac3D; 
end
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); p.nt=pde.grid.nElements; 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t; p.mesh.be=e;  % background mesh  
p.fuha.e2rs=@e2rs; % and refinement selector 
p.u=zeros(p.np,1); p.u=[p.u; par]; p=setfemops(p); p.sw.para=2; 
p.usrlam=1:1:4; p.nc.lammax=5;  p.sol.ds=0.2; p.nc.dsmax=0.2; p.nc.nsteps=20; 
p.nc.sf=1e3; % stiff spring constant for DBC via Robin-BC 