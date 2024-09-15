function p=shinit(p,lx,ly,nx,ny,par,Lfn) % init with prep. of mesh-adaption and usrlam
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.Lfn=Lfn; % filename for matrices 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
p.np=nx*ny; p.lx=lx; p.ly=ly; p.nx=nx; p.ny=ny; p.nu=p.np; p.sol.xi=1/(p.nu);
p.u=zeros(p.np,1); p.u=[p.u; par']; % initial guess, and parameters 
p.usrlam=-0.5:0.5:1; % compute point and write to disk at these par values
p=oosetfemops(p);    % generate diff matrices 
p.nc.nsteps=20; p.sw.foldcheck=0; p.nc.neig=4; p.nc.eigref=-0.1; 
p.plot.auxdict={'\lambda','c2','c3'}; p.plot.pstyle=-1; 
p.nc.ilam=1; p.nc.lammin=-1; p.nc.lammax=1; 
p.sol.ds=0.1; p.nc.dsmax=0.1; p.nc.dsmin=0.01; 