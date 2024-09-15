function p=shinit(p,lx,nx,par) % SH-init, matrix-free dct version 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.plot.pstyle=-1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
p.np=nx; p.lx=lx; p.nu=p.np; p.sol.xi=1/(p.nu);
p.u=zeros(p.np,1); p.u=[p.u; par']; %del=1e-6; p.usrlam=(-0.5:0.5:1)+del; 
p=oosetfemops(p);    % generate diff matrices 
p.nc.nsteps=20; p.sw.foldcheck=0; p.nc.foldtol=0.01; 
p.plot.auxdict={'\lambda','c2','c3'}; 
p.nc.ilam=1; p.nc.lammin=-1; p.nc.lammax=2; p.sol.ds=0.1; p.nc.dsmax=0.1; 
p.sw.bifcheck=2; p.nc.mu1=0.5; p.nc.mu2=0.02; p.sw.verb=2; 