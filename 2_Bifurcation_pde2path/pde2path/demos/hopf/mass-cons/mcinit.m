function p=mcinit(p,lx,nx,par,ndim) 
p=stanparam(p); screenlayout(p); % set standard parameters and screenlayout 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.sw.jac=0; % rhs and Jac
p.nc.neq=2; p.nc.ilam=5; p.fuha.outfu=@hobra; % number of eq, cont-param, output
p.sw.verb=2; p.sw.spcalc=1; 
switch ndim % domain and BC, depending on spatial dim 
  case 1; pde=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; p.x0i=10; % index for ts-plot
        bc=pde.grid.neumannBC('0'); p.plot.pcmp=[1 2];% BC 
  case 2; pde=stanpdeo2D(lx,lx/2,2*lx/nx); p.vol=2*lx^2; p.x0i=30; 
        bc=pde.grid.neumannBC('0'); 
end 
p.pdeo=pde; p.sw.sfem=-1; p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;  
p=setfemops(p); % setfemops calls oosetfemops in problem dir 
al=par(5); beta=par(6); u0=-beta/2-sqrt(beta^2/4+al); %u0=-1*sqrt(al+beta); 
u=u0*ones(p.np,1); v=-u; p.u=[u;v; par]; % initial guess and pars 
p.usrlam=[];  % user-vals for output
p.plot.cm=cool; p.plot.bpcmp=10; % colormap for soln plot, component for branch-plot 
[p.u,res]=nloop(p,p.u);fprintf('first res=%g\n',res); % start-point for cont 
p.file.smod=10; p.sol.ds=0.1; p.nc.dsmax=0.5; % saving, stepsize, max stepsize 
p.sw.bifcheck=2; p.nc.neig=20; % method for bifcheck, and # Evals used 
p.nc.mu1=1; p.nc.mu1=0.5; % be relaxed about possible bif-detection
p.sol.ds=-0.1; p.nc.dsmax=1; p.nc.dlammax=1; p.nc.tol=1e-6; 