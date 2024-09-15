function p=acinit(lx,nx,par,femorder) % init with prep. of mesh-adaption and usrlam
p=stanparam; screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.e2rs=@e2rs; 
pde=stanpdeo1D(lx,2*lx/nx); p.pdeo=pde; % standard domain and mesh
ne=pde.grid.nPoints-1; 
p.hofem.femorder=femorder; p=two2lob(p); % convert to lobatto mesh 
p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu);
p.u=zeros(p.np,1); p.u=[p.u; par']; % initial guess, and parameters 
p.usrlam=-0.5:0.5:1; % compute point and write to disk at these par values
p=setfemops(p);    % generate FEM matrices 
p.nc.nsteps=20; p.sw.foldcheck=1; p.mesh.maxt=100; 
p.plot.auxdict={'c','\lambda','gamma','d'}; 
p.nc.ilam=2; p.nc.lammin=-0.2; p.nc.lammax=1; 
p.sol.ds=0.1; p.nc.dsmax=0.1; 