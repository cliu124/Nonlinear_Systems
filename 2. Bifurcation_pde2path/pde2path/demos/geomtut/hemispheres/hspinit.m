function p=hspinit(nx,par) % hemisphere init 
p=stanparam(); p=stanparamX(p);  % set stanparam, adapt to X; next reset some: 
p.fuha.sG=@sGhs; p.fuha.sGjac=@sGjac; p.sw.spcalc=0; p.sw.bifcheck=0; 
%pde=diskpdeo2(1,nx,round(nx/2)); % disk preimage discretization, pde-object 
pde=diskpdeo2(1,nx,round(nx/2));
% not stored, only p.DBC, p.tri and p.X (generated below) used subsequently
p.np=pde.grid.nPoints; p.nu=p.np; p.nt=pde.grid.nElements;  % store dimensions
p.sol.xi=1/p.nu; p.n0=p.np; % u-vs-lam weight, initial mesh size (for coarsening) 
p.nc.neq=1; p.sw.jac=0; % here, for simplicity, numerical Jacs 
po=pde.grid.p; x=po(1,:)'; y=po(2,:)'; u=0*ones(p.np,1); % set ICs 
p.u=[u; par]; p.X=[x,y,sqrt(1+1e-12-x.^2-y.^2)]; 
p.tri=pde.grid.t(1:3,:)'; % store initial X and tri 
p.NEU=pde.grid.e(1:2,:)'; p.idx=unique(p.NEU(:)); p.X(p.idx,3)=0; % edges/points for BCs 
p=oosetfemops(p); p.plot.auxdict={'H','V','A'}; 
p.u(p.nu+2)=getV(p,p.u); p.u(p.nu+3)=getA(p,p.u); % get initial vol & area
p.mpos=1; % use M for computing center of gravity 