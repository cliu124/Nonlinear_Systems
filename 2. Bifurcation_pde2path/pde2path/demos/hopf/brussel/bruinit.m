function p=bruinit(p,lx,nx,par,ndim) % init for Brusselator 
p=stanparam(p); screenlayout(p); p.nc.ilam=2; p.nc.neq=3; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.outfu=@hobra; 
p.fuha.e2rs=@e2rsbru; p.nc.sig=0.1; p.nc.ngen=2; % meshada settings 
switch ndim
    case 1;pde=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; 
        bc=pde.grid.neumannBC('0'); p.x0i=10; % index for ts-plot
    case 2; ly=lx/4; %ly=lx; 
        pde=stanpdeo2D(lx,ly,2*lx/nx); p.vol=4*lx*ly; 
        bc=pde.grid.neumannBC('0'); p.x0i=30; p.plot.pstyle=2; 
end 
pde.grid.makeBoundaryMatrix(bc); p.nc.sf=1e3; 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu; 
a=par(1); b=par(2); c=par(3); d=par(4); % set Iguess depending on params 
u=a*ones(p.np,1); v=b*u/a^2; w=c*u/d; p.u=[u;v;w]; p.u(p.nu+1:p.nu+7)=par; 
p.sol.ds=0.025; p.nc.dsmax=0.025; p.nc.dsmin=1e-8; p.file.smod=5; p.plot.cm=hot;
p.nc.imax=10; p.sw.para=1; p.nc.dsinciter=p.nc.imax/2; p.plot.bpcmp=9; 
p=setfemops(p); 