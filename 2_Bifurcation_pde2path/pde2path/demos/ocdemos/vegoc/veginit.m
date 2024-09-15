function p=veginit(p,lx,ly,nx,sw,rho,ndim) % init-routine for vegOC 
p=stanparam(p); p.nc.neq=4; p.fuha.sG=@vegsG; p.fuha.jcf=@vegjcf; 
p.fuha.outfu=@ocbra; p.ndim=ndim; p.d1=0.05; p.d2=10; h=2*lx/nx; 
if ndim==1; pde=stanpdeo1D(lx,h); p.vol=2*lx; p.plot.pstyle=5; 
else; s.sym=2; 
    pde=stanpdeo2D(lx,ly,h,s); p.vol=4*lx*ly; p.plot.pstyle=2; 
end 
p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; 
p.sol.xi=0.005/p.np; p.pdeo=pde; p.sw.sfem=-1; p=setfemops(p); 
p.xplot=lx; p.sw.spcalc=0; p.sw.jac=1; p.file.smod=100; 
p.nc.dsmin=1e-6; p.nc.lammin=1e-10; p.nc.nsteps=150; p.plot.bpcmp=0; 
p.nc.ilam=8; p.nc.dlammax=5; p.sw.jac=0; p.nc.imax=60;
switch sw % build initial guess 
    case 1; % rho=0.03  
    p.usrlam=[4 10 20 26 28]; 
    par=[rho 1e-3 0.5 0.03 0.005 0.9 1e-3 34 0.01 0.1 1 1.1 0.3];  
    foss=[400 12 0.5 0.9]; p.sol.ds=-0.5; u=foss'*ones(p.np,1)'; 
    u0=[u(1,:) u(2,:) u(3,:) u(4,:)]; p.u=u0(:);   p.nc.dsmax=1; 
    p.nc.lammax=34;
    case 2; % rho=0.03, larger R (for comp. with uncontrolled case) 
    p.usrlam=[10 20 30 40 50 60 70 80 90 100 110 120 130]; 
    par=[rho 1e-3 0.5 0.03 0.005 0.9 1e-3 130 0.01 0.1 1 1.1 0.3];  
    foss=[1500 50 0.5 0.9]; p.sol.ds=-0.5; u=foss'*ones(p.np,1)'; 
    u0=[u(1,:) u(2,:) u(3,:) u(4,:)]; p.u=u0(:); p.nc.dsmax=10; 
    p.nc.lammax=130;
    case 3; % rho=0.03  (more usr-lam-values) 
    p.usrlam=[4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36]; 
    par=[rho 1e-3 0.5 0.03 0.005 0.9 1e-3 24 0.01 0.1 1 1.1 0.3];  
    foss=[400 12 0.5 0.9]; p.sol.ds=-0.5; u=foss'*ones(p.np,1)'; 
    u0=[u(1,:) u(2,:) u(3,:) u(4,:)]; p.u=u0(:);   p.nc.dsmax=1; 
end
p.u=[p.u; par']; r=resi(p,p.u); fprintf('inires=%g\n',norm(r,'inf')); 
[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res); 
plotsol(p,1,1,p.plot.pstyle); p.nc.imax=10;
