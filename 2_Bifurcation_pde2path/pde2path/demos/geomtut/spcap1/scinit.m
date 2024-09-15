function p=scinit(nx,par) % spherical cap, init 
p=stanparam(); p=stanparamX(p);  % set stanparam, adapt to X; then reset some
p.fuha.sG=@sGsc; p.fuha.sGjac=@scjac; p.sw.spcalc=0; p.sw.bifcheck=0; 
pde=diskpdeo2(1,nx,round(nx/2)); % disk preimage discretization, pde-object 
% not stored, only p.DBC, p.tri and p.X (generated below) used subsequently
p.np=pde.grid.nPoints; p.nu=p.np; p.nt=pde.grid.nElements;  % store dimensions
p.sol.xi=1/p.nu; p.n0=p.np; % u-vs-lam weight, initial mesh size (for coarsening) 
p.nc.neq=1; p.sw.jac=0; % here, for simplicity, numerical Jacs 
po=pde.grid.p; x=po(1,:)'; y=po(2,:)'; u=0*ones(p.np,1); % set ICs 
p.u=[u; par]; p.X=[x,y,0*x]; p.tri=pde.grid.t(1:3,:)'; % store initial X and tri 
p.DIR=pde.grid.e(1:2,:)'; p.idx=unique(p.DIR(:)); % edges, and points for BCs 
p=oosetfemops(p); p.plot.auxdict={'H','V','A'}; % dummy oosetfemops! 
p.u(p.nu+2)=getV(p,p.u); p.u(p.nu+3)=getA(p,p.u); % get initial vol & area
p.nc.lammax=200; p.nc.dsmax=3; p.nc.dlammax=3; p.file.smod=10;