function p=shinit(p,lx,ly,nx,ny,par) % init with prep. of mesh-adaption and usrlam
p=stanparam(p); screenlayout(p); p.nc.neq=1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
p.np=nx*ny; p.lx=lx; p.ly=ly; p.nx=nx; p.ny=ny; p.nu=p.np; p.sol.xi=1/(p.nu);
p.u=zeros(p.np,1); p.u=[p.u; par']; % initial guess, and parameters 
p=oosetfemops(p);    % generate diff matrices 
p.nc.nsteps=20; p.sw.foldcheck=0; p.nc.neig=4; p.nc.eigref=0; 
p.plot.auxdict={'\lambda','c2','c3'}; p.plot.pstyle=-1; 
p.nc.ilam=1; p.nc.lammin=-1; p.nc.lammax=1; 
p.sol.ds=0.1; p.nc.dsmax=0.1; p.nc.dsmin=0.01; 
p.nc.mu1=1; p.nc.mu2=0.05; % allow rather crude Eval detection and loc. 