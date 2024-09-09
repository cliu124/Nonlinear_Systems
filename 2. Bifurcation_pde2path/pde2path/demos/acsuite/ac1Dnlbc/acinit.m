function p=acinit(p,lx,nx,par) 
p=stanparam(p); screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; % domain and mesh
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
[po,t,e]=getpte(p); p.mesh.bp=po; p.mesh.bt=t;  
p.mesh.be=e; p.mesh.nt=size(t,2); % background-mesh (for oomeshadac)
p.u=zeros(p.np,1); p.u=[p.u; par']; % initial guess, and parameters 
p=setfemops(p);    % generate FEM matrices 
p.nc.nsteps=20; p.sw.foldcheck=1; p.mesh.maxt=100; 
p.plot.auxdict={'c','lambda','gamma','d','alpha','beta'}; 