function p=acinit(p,nx) % init-routine 
p=stanparam(p); p.nc.neq=1; p.nc.ilam=4; 
p.fuha.G=@acf; p.fuha.Gjac=@acjac; 

% set geometry to rectangle and BC to hom. Dirichlet (stiff spring) 
lx=1; ly=0.9; [p.mesh.geo, bc]=recdbc1(lx,ly,1e3);
p.fuha.bc=@(p,u) bc; p.fuha.bcjac=@(p,u) bc;
% generate standard mesh and init this as the base mesh for refinement
ny=nx; p=stanmesh(p,nx,ny); p=setbmesh(p); 

% initial condition, stepping and parameter settings
p.u=0*ones(p.nu,1); p.sol.ds=0.5; p.nc.nsteps=10; p.nc.lammax=10;
p.sol.xi=1/p.np; p.nc.amod=0; p.nc.lamdtol=0.8; 

% initialize auxiliary variables, here parameters of PDE
par(4)=1;    % linear coefficient of f
par(1)=0.25; par(2)=-0.2; par(3)=0.01; % diffusion const, lin and nonlin 
p.u=[p.u ;par']; % append pars to u 