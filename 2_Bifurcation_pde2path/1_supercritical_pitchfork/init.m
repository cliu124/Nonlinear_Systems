function p=init(p,par)  % init routine for AC on interval with pBC 
p.np=2; %dimension of ODE systems. 

%% default setup, no need to change. 
p=stanparam(p); screenlayout(p); p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
pde=stanpdeo1D(1,4/p.np); p.pdeo=pde; % domain and mesh
p.nu=p.np; p.sol.xi=1/(p.nu); [po,t,e]=getpte(p);
p.mat.M=eye(p.np);

%% define the trivial solution. In most cases just zero. 
p.u=zeros(p.np,1); 
p.u=[p.u; par']; % initial guess (here 0, explicitly known) 

p.sw.foldcheck=1; 
p.plot.auxdict={'lambda'}; %for plotting to show lambda
p.plot.pstyle=1; %defaulty
p.nc.nsteps=200; %default continuation step
p.sw.bifcheck=2; %bifcheck=2: check the bifurcation points. 
p.nc.ilam=1; %select ilam th parameter as the continuation parameter. For only one parameter system, this is just 1. 
p.nc.lammax=2; %the maximum lambda you want to continue to
p.sol.ds=0.001; %the initial step size of numerical continuation
p.nc.dsmax=0.05; %the maximal step size of numerical continuation
