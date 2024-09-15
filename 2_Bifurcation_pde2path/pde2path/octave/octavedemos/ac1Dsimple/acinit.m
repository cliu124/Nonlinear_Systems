
function p=acinit(p,lx,nx,par) % init AC as a rather generic scalar PDE 
p=stanparam(p); screenlayout(p); % standard settings/screenlayout 
p.nc.neq=1; p.sw.sfem=-1; % scalar problem, implementation with OOPDE: 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; % function handles to rhs and Jacobian
pde=stanpdeo1D(lx,2*lx/nx); % standard PDE object 1D, yields domain and mesh
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
p.u=zeros(p.np,1); p.u=[p.u; par']; % initial sol, with param appended 
p=oosetfemops(p);    % generate FEM matrices 
p.nc.nsteps=20; % # of continuation steps 
p.nc.ilam=2; p.nc.lammax=2; % continue in par(2)=lambda, up to lammax 
p.sol.ds=0.1; p.nc.dsmax=0.1; % starting and max value for contin. step size
p.sw.foldcheck=1; % detect and localize folds 
p.plot.auxdict={'c','lambda','gamma'}; % parameter names 
p.plot.fancybd=0; 