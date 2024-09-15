function p=slinit(p,lx,ly,nx,sw,ndim) % init-routine 
p=stanparam(p); % set generic parameters to standard, if needed reset below..
p.fuha.sG=@slsG; p.fuha.sGjac=@slsGjac; p.fuha.outfu=@ocbra; % rhs and branch
p.fuha.jcf=@sljcf; p.fuha.con=@slcon; % current-val, and fun to get k from u

h=2*lx/nx;  % mesh-size; now generate domain and mesh depending on dimension 
if ndim==1; pde=stanpdeo1D(lx,h); p.vol=2*lx; % 1D
else pde=stanpdeo2D(lx,ly,h); p.vol=4*lx*ly;  % 2D 
end 

p.pdeo=pde; p.np=pde.grid.nPoints; % save the pde object (and some more) in p 
p.sol.xi=1/p.np; p.nc.neq=2; p.nu=p.np*p.nc.neq;  % weight, #eqns, and DoF 
p.sw.sfem=-1; p=setfemops(p); % use OOPDE, generate FEM matrices 
p.sw.spcalc=0; p.file.smod=100; % don't compute EVals, save each 100th point 
par=[0.03; 0.55; 0.5; 0.5]; p.nc.ilam=2; p.sol.ds=0.1;   % startup param 
% r=par(1); bp=par(2); cp=par(3); D=par(4);  % to recall the meaning of par
p.usrlam=[0.55 0.6 0.65 0.7 0.75]; % target-values for lam 
p.nc.dsmin=1e-6; p.nc.dsmax=0.5; p.nc.lammax=0.8; p.nc.lammin=0.549; 

switch sw  % choose initial guess according to switch 
  case 1; u=0.3*ones(p.np,1); v=-13*ones(p.np,1); u0=[u v]; p.u=u0(:); % FSC
  case 2; u=2*ones(p.np,1); v=-4*ones(p.np,1); u0=[u v]; p.u=u0(:); % FSM
end

p.u=[p.u; par]; % append pars, find 1st sol from initial guess then plot
[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res); plotsol(p,1,1,1); 