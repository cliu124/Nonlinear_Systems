function p=ac6init(p) 
% init-routine 
p=stanparam(p); % set generic parameters to standard, if needed reset below..
dir=sprintf('%s',inputname(1)); p=setfn(p,dir); screenlayout(p); 
% PDE standard implementation
p.nc.neq=1; p.fuha.G=@ac_G; p.fuha.Gjac=@ac_Gjac;  
p.fuha.bc=@ac6bcfx; p.fuha.bcjac=@ac6bcfx_jac; 

lx=0.5;ly=0.5; p.mesh.geo=hexgeo(lx,ly); p=stanmesh(p,0.1); p=setbmesh(p); % geometrie and mesh
p.nc.nsteps=50; p.sol.ds=0.05; p.sol.xi=1/p.np; p.nc.lammin=-1; p.nc.lammax=2; 

p.u=0*ones(p.nu,1);    % initial point 
% parameters
par(1)=0;    % linear cofficient of f
par(2)=0.25; % diffusion coefficient
par(3)=1;    % nonlinear coefficient of f
p.u = [p.u; par']; p.nc.ilam=1;

