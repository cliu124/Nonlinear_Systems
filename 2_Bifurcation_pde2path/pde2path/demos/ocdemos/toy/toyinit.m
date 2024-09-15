function p=toyinit(p,nx,nt,par) % init file for the CPS
lx=1; p=stanparam(p); screenlayout(p); % set standard parameters  
p=hostanparam(p); 
p.hopf.tom.Stats='on'; p.hopf.tom.Stats_step='on'; p.hopf.tom.Monitor=3; p.hopf.tom.Order=2;
p.hopf.tom.AbsTol=1e-04; p.hopf.tom.Itnlmax=20; p.hopf.tom.Itlinmax=20;
p.hopf.tom.FJacobian=@hojac; p.hopf.tom.BCJacobian=@hobcjac; p.hopf.tom.lu=0; p.hopf.tom.vsw=0; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.sw.jac=1; % rhs and Jac
p.nc.neq=4; p.nc.ilam=1;
pde=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; p.x0i=10; % index for ts-plot
p.nc.sf=1e3; p.pdeo=pde; % OOPDE setting of BC 
p.sw.sfem=-1; p.np=1; p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu;  
p.sol.ds=0.1; p.nc.dsmax=0.5; % saving, stepsize, max stepsize
p.hopf.T=2*pi; p.hopf.tom.Nmax=1e10; p.hopf.lam=1;
p.mat.M=eye(4,4); % replace mass and "laplacian"-Matrix
p.mat.K=0*eye(4,4); p.hopf.tom.M=p.mat.M;
p.u=[zeros((nx-1)*4,1);par]; % Dummy solution, only paras are necessary
p.hopf.tl=nt; p.hopf.tau=0; p.hopf.muv1=0; p.hopf.muv2=0;
p.hopf.y0dsw=2; p.sw.jac=1; p.hopf.jac=1;
% construct CPS
p.hopf.t=linspace(0,1,nt);
p.hopf.y=[cos(p.hopf.T*p.hopf.t);sin(p.hopf.T*p.hopf.t);ones(1,nt);zeros(1,nt)];
p.hopf.y0d=1/(2*pi)*[(-2*pi*sin(p.hopf.T*p.hopf.t));2*pi*cos(p.hopf.T*p.hopf.t);zeros(1,nt);zeros(1,nt)];
y=p.hopf.y; yshift=circshift(y(:,1:end-1),-1,2); p.hopf.y=[yshift yshift(:,1)];
p.hopf.T=2*pi/par(3);