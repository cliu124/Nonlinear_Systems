function p=secoinit(lx,nx,par) % init for MWBS'97 semiconductor model 
p=stanparam(); p.nc.neq=2;     % init with stanparam, 2-compo-system 
p.sol.ds=-0.01; p.sol.dsmax=0.05; p.sw.bifcheck=2; % reset some pars 
p.fuha.outfu=@secobra; p.plot.bpcmp=8; % output function handle, and compo 
p=setbel(p,0,1e-6,10,@lss); % use bordered elim. lss (always good in 1D), 
% with border width 0, tolerance 1e-6, and at most 10 iterations 
pde=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; % standard 1D PDE object 
p.np=pde.grid.nPoints; p.pdeo=pde; % store number of grid-points, and pdeo 
p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu; % DoFs, and weight for arclength 
p.nc.lammin=0; p.nc.lammax=5; p.nc.dsmax=0.1; % lam-range, and max stepsize
j0=par(1); tau=par(3); als=j0/(tau*(j0^2+1)); us=als+j0; % initial sol
u=als*ones(p.np,1);v=us*ones(p.np,1);p.u=[u;v;par]; % append pars and store 
p.nc.ilam=1; % select the (initial) primary active par, here j0 
p.sw.sfem=-1; p=oosetfemops(p); % set FEM matrices 