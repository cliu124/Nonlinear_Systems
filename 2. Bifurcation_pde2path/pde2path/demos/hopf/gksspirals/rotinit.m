function p=rotinit(p,nx,par) % init for gksspirals-demo, legacy sfem=1 setting
p=stanparam(p); screenlayout(p); 
p.file.dircheck=0; p.nc.ilam=1; p.nc.neq=2; p.fuha.outfu=@hobra; 
p.fuha.sG=@sG;p.fuha.sGjac=@sGjac; % rhs
p.eqn.c=isoc([[0.01,0];[0,0.015]],2,1); % diffusion tensor (for setfemops) 
p.eqn.b=0; p.eqn.a=0; % no conv. or linear terms (put these into nodalf) 
p.fuha.bc=@nbc; p.fuha.bcjac=@nbc; % BCs 
p.mesh.geo=circgeo(1,nx); hmax=2/nx; p=stanmesh(p,hmax); 
p.sw.sfem=1; p.vol=2*pi; p.sw.bifcheck=2; p.nc.neig=20; 
p=setfemops(p); % here using legacy setfemops, diff
p.nu=p.np*p.nc.neq; p.sol.xi=1/p.nu; p.plot.pstyle=3; 
u=0*ones(p.np,1); v=u; p.u=[u;v]; p.u(p.nu+1:p.nu+length(par))=par; 
p.sol.ds=0.1; p.nc.dsmin=1e-8; p.nc.dsmax=0.1; p.file.smod=5; p.plot.cm=hot;
p.nc.imax=10; p.sw.para=1; p.nc.dsinciter=p.nc.imax/2; p.plot.bpcmp=0; 
r=resi(p,p.u); fprintf('inires=%g\n',norm(r,p.sw.norm));
[p.u,res]=nloop(p,p.u);fprintf('first res=%g\n',res); p.sw.jac=1; 