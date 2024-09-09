function p=acinit(p,lx,nr,na,par) % init with prep. of mesh-adaption and usrlam
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.plot.pstyle=-1; 
p.sw.sfem=-1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.nr=nr; p.na=na; p.lx=lx; 
p=oosetfemops(p); p.Om=pi*p.lx^2*(p.r(1)^2-p.r(end)^2); 
p.nu=p.np; p.sol.xi=1/(p.nu);
p.u=zeros(p.np,1); p.u=[p.u; par']; % initial guess, and parameters 
%plotsol(p); pause 
p.usrlam=-0.5:0.5:1; % compute point and write to disk at these par values
%p=oosetfemops(p);    % generate FEM matrices 
p.nc.nsteps=20; p.sw.foldcheck=1; p.mesh.maxt=100; 
p.plot.auxdict={'c','\lambda','gamma','d'}; 
p.nc.ilam=2; p.nc.lammin=-2; p.nc.lammax=20; 
p.sol.ds=0.1; p.nc.dsmax=0.1; 