function p=cGLinit(p,lx,nx,par,ndim) % (generic) init routine for cGL problem 
p=stanparam(p); screenlayout(p); % set standard parameters and screenlayout 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.sw.jac=1; % rhs and Jac
p.nc.neq=2; p.nc.ilam=1; p.fuha.outfu=@hobra; % number of eq, cont-param, output
p.usrlam=[-0.25 -0.2 -0.1 0 0.5  1 2 3];  % user-vals for output
p.plot.cm=hot; p.plot.bpcmp=9; % colormap for soln plot, comp.for branch-plot 
switch ndim % domain and BC, depending on spatial dim 
  case 1; pde=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; p.x0i=10; % index for ts-plot
        bc=pde.grid.neumannBC('0'); % BC 
  case 2; pde=stanpdeo2D(lx,lx/2,2*lx/nx); p.vol=2*lx^2; p.x0i=30; 
        bc=pde.grid.robinBC('1','0'); 
  case 3; pde=stanpdeo3D(lx,lx/2,lx/4,2*lx/nx); p.vol=0.5*lx^3; p.x0i=200; 
        bc=pde.grid.robinBC('1','0'); 
        p.plot.ng=20; p.plot.lev=[-0.1 0.1]; % 3D specific plot settings 
        p.plot.levc={'blue','red'}; p.plot.alpha=0.5;    
  case 4; sw.sym=2; pde=stanpdeo2D(lx,lx,nx,nx,sw); p.vol=4*lx^2; p.x0i=30; 
        bc=pde.grid.neumannBC('0'); p.nc.lammax=2; p.usrlam=[1 2]; p.plot.cm='cool'; 
        p.hopf.ax='unif'; p.hopf.lay=[1 4]; p.hopf.pind=1:5:16; 
        p.hopf.xtics=[]; p.hopf.ytics=[]; 
end 
pde.grid.makeBoundaryMatrix(bc); p.nc.sf=1e3; p.pdeo=pde; % OOPDE setting of BC 
p.sw.sfem=-1; p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;  
p=setfemops(p); % setfemops calls oosetfemops in problem dir 
u=0*ones(p.np,1); v=u; p.u=[u;v; par]; % initial guess (here trivial) and pars 
[p.u,res]=nloop(p,p.u);fprintf('first res=%g\n',res); % start-point for cont 
p.file.smod=10; p.sol.ds=0.1; p.nc.dsmax=0.5; % saving, stepsize, max stepsize 
p.sw.bifcheck=2; p.nc.neig=20; % method for bifcheck, and # Evals used 
p.nc.mu1=0.5; % be relaxed about possible bif-detection