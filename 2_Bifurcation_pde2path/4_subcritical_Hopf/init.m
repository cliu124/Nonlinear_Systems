function p=init(p,par)  % init routine for AC on interval with pBC 
p.np=2; %dimension of ODE systems. 

p=stanparam(p); screenlayout(p); p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
pde=stanpdeo1D(1,4/p.np); p.pdeo=pde; % domain and mesh
 p.nu=p.np; p.sol.xi=1/(p.nu); [po,t,e]=getpte(p);
p.mat.M=eye(p.np);
p.u=zeros(p.np,1); p.u=[p.u; par']; % initial guess (here 0, explicitly known) 
p.nc.nsteps=20; p.sw.foldcheck=1; p.plot.auxdict={'lambda'}; 
p.plot.pstyle=1; p.nc.nsteps=100; %p.sw.jac=1; 
p.sw.bifcheck=2; p.nc.ilam=1; p.nc.lammax=2; p.sol.ds=0.001; p.nc.dsmax=0.1;
