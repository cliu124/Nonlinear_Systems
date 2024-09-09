function p=acfront_init(p) 
% initialization routine for "acfront" demo
% see manual or routine "stanparam" for explanations of settings

% standard settings, window arrangement, file names
p=stanparam(p); dir=sprintf('%s',inputname(1)); p=setfn(p,dir);

% PDE standard implementation
p.nc.neq=1;           % set number of equations for pde
p.fuha.G=@acfront_G;         % set functions f in pde
p.fuha.Gjac=@acfront_Gjac;  % set jacobian of f (not required if p.jsw=0)

% domain and mesh
lx=0.1; ly=2; [p.mesh.geo,bc]=recnbc1(lx,ly); 
p.fuha.bc=@(p,u) bc; p.fuha.bcjac=@(p,u) bc;
nx=2; ny=50; p=stanmesh(p,nx,ny); p=setbmesh(p); 

% continuation settings
p.nc.ilam=1;       % set primary parameter index (here no other active parameters)
p.sol.ds=0.5; p.nc.nsteps=30; p.nc.lammax=20;
p.sol.xi=1/p.np; p.nc.lamdtol=0.8; 

% initial condition, stepping and parameter settings
p.u=0*ones(p.nu,1); 
par(1)=0.1; % linear cofficient of f
par(2)=1; % second nonlinear (energy well diff.)
par(3)=0; % advection term
p.u = [p.u; par'];
p=setfemops(p); % generate FEM operators for this


