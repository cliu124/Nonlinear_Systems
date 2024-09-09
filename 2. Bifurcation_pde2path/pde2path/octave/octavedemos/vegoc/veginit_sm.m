function p=veginit_sm(p,lx,ly,nx,sw,rho,ndim) % init-routine for vegOC 
p=stanparam(p); p.fuha.jcf=@vegjcf; p.fuha.outfu=@ocbra; p.ndim=ndim; 
p.nc.neq=4; p.fuha.sG=@vegsG; p.sw.jac=0; % #eqns, rhs, numerical jacs 
p.d1=0.05; p.d2=10; h=2*lx/nx; p.sw.sfem=-1; % diff. constants, mesh-width
if ndim==1; pde=stanpdeo1D(lx,h); p.vol=2*lx; p.plot.pstyle=5; % 1D pde-object
else pde=stanpdeo2D(lx,ly,h); p.vol=4*lx*ly; p.plot.pstyle=2;  % 2D
end 
p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; p.sol.xi=0.005/p.np; 
p.pdeo=pde; p=setfemops(p); % create FEM matrices 
p.sw.spcalc=0; p.file.smod=100; p.nc.ilam=8; p.nc.dlammax=5; % some settings
p.nc.dsmin=1e-6; p.nc.lammin=1e-10; p.nc.nsteps=150; p.plot.bpcmp=1; 
switch sw % build initial guess 
    case 1; % rho=0.03  
    p.usrlam=[4 10 20 26 28]; % desired user values for R 
    par=[rho 1e-3 0.5 0.03 0.005 0.9 1e-3 34 0.01 0.1  1  1.1 0.3];  % param., i.e., 
    %    rho   g  eta  d    del beta  xi   R  r_u r_w  c  p   alpha
    fcss=[400 12 0.5 0.9]; % guess for (starting) FCSS (flat canonical steady state)  
    p.sol.ds=-0.5; u=fcss'*ones(p.np,1)'; u0=[u(1,:) u(2,:) u(3,:) u(4,:)]; 
    p.u=u0(:); p.nc.lammax=34; p.nc.dsmax=1; 
   % other cases omitted here
end
p.u=[p.u; par']; % append params to FCSS guess and Newton to true soln
r=resi(p,p.u); fprintf('inires=%g\n',norm(r,'inf')); 
[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res); plotsol(p,1,1,p.plot.pstyle); 