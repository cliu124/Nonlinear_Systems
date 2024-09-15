function p=bratu_init(p,lx,ly,nx) 
% initialization routine for "bratu" demo
% see manual or routine "stanparam" for explanations of settings
% standard settings, window arrangement, file names
p=stanparam(p); dir=sprintf('%s',inputname(1)); p=setfn(p,dir); screenlayout(p);

% PDE standard implementation
p.nc.neq=1;        % set number of equations for pde
p.fuha.G=@bratu_G;       % set functions f in pde
p.fuha.Gjac=@bratu_Gjac;   % set jacobian of f (not required if p.jsw=3)
p.nc.nq=0;         % number of auxiliary equations

% domain and mesh
[p.mesh.geo, bc]=recnbc1(lx,ly); p.fuha.bc=@(p,u) bc; p.fuha.bcjac=@(p,u) bc;
ny=round(nx*ly/lx);p=stanmesh(p,nx,ny); p=setbmesh(p); % p=setfemops(p);

% continuation settings
p.nc.ilam=1;       % set primary parameter index (here no other active parameters)
p.sol.xi=1/p.np; p.nc.nsteps=50; p.nc.dsmax=1; p.nc.dlammax=0.05; 
p.nc.lammin=0.05; p.sol.ds=0.05;

% initial point and parameters
p.u=0.1*ones(p.np,1);
par(1)=0.2;  % parameter 1 
par(2)=1;    % parameter 2
p.u =[p.u; par'];

% newton step to get good initial point
[p.u,p.sol.res,p.sol.iter]=nloop(p,p.u); % correct initial guess 
fprintf('first residual=%g with %i iterations\n',p.sol.res,p.sol.iter); % and inform user 
plotsol(p,1,1,p.plot.pstyle); % plot this guess

