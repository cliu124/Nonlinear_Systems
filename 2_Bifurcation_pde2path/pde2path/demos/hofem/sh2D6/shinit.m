function p=shinit(p,lx,ly,nx,ny,dsw,par,nref,hofem) % Swift-Hohenberg as 2 component system, with 
% singular M=diag(M,0) mass-matrix to compute correct eigenvalues 
p=stanparam(p); p.nc.neq=2; p.ndim=2; p.nc.neig=4; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
p.fuha.outfu=@shbra; p.sw.spcalc=0; p.plot.bpcmp=0; 
p.sol.ds=0.01; p.sol.dsmax=0.1; p.pm.resfac=1e-3; p.sw.bifcheck=2; 
p.sw.spcont=0; p.sw.spcalc=1; p.sw.sfem=-1; p.sw.foldcheck=1; 
p.hofem=hofem; 
switch dsw
 case 1; % create stanpdeo2D, convert to 6node triangulation, 
   sw.sym=1; pde=stanpdeo2D(lx,ly,nx,ny, sw); p.pdeo=pde; p.hofem.t2sinterpol=0; p=tri2six(p); 
 case 2;  % directly create 6node triangulation 
   pde=stanpdeo2D(lx,ly,2,2); % dummy pdeo, filled below 
   [p.nt,p.np,p.po,p.tri,p.efl,p.gfl]=trgl6_rec(lx,ly,nx,nref); % 6-node-triangulation 
   c3=c3fu(p.tri,p.nt); size(c3), c3=[c3, ones(size(c3,1),1)]; % convert to 3-node triangle-mesh 
   pde.grid.p=p.po'; pde.grid.t=c3'; p.pdeo=pde; 
end
p.plot.pstyle=2; p.plot.cm='cool'; p.plot.axis='image'; 
p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;
p.nc.lammin=-4; p.nc.lammax=2; p.nc.ds=0.01; p.nc.dsmax=0.1; 
u=0*ones(p.np,1); v=u; u0=[u v]; p.u=u0(:); 
p.u=[p.u; par']; p.nc.ilam=1; p.file.smod=5; 
p=setfemops(p); p.nc.nsteps=100; 
p.Om=sum(sum(p.mat.M(1:p.np,1:p.np),1)); 