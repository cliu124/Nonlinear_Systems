function p=cheminit(p,lx,ly,nx,ny) 
% init-routine 
p=stanparam(p); p.dim=2; p.sw.sfem=-1; huclean(p); 
p.nc.neq=2; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.outfu=@chembra; 
p.fuha.cfu=@cfu; 
sw.sym=2; pde=stanpdeo2D(lx,ly,nx,ny,sw);
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=2*p.np; p.nt=pde.grid.nElements; 
p.sol.xi=1/p.nu; p.nc.dsmin=1e-4; p.nc.dsmax=1; 
p.nc.dlammax=0.5; p.nc.neig=min(20,p.np); 
p.sw.spcalc=0; p.plot.pstyle=2; p.plot.pcmp=1;
p.plot.pcmp=2; p.sw.para=1; p.sol.ds=0.1;
p.nc.dlammax=2; p.nc.lammax=25; 
u=1*ones(p.np,1); v=0.5*ones(p.np,1); p.u=[u; v; 10]; p.nc.ilam=1; 
p=oosetfemops(p); plotsol(p,1,1,1); p.vol=4*lx*ly; % needed in chembra 
