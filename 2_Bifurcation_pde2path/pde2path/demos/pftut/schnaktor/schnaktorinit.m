function p=schnaktorinit(p,dom,nx,par,psw,varargin)
p=stanparam(p); screenlayout(p); p.nc.neq=2; p.sw.sfem=-1;
p.plot.auxdict={'\lambda','\sigma','d','||u_1||_{\infty}','min(|u_1|)'};
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.spjac=@spjac; 
p.fuha.qf=@qf; p.fuha.qfder=@qfder;
ny=round(dom(2)/dom(1)*nx); lx=dom(1); ly=dom(2); p.vol=4*lx*ly; 
pde=stanpdeo2D(lx,ly,nx,ny,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; p.nc.neig=30; p.nc.nsteps=50; 
p.sol.xi=1/p.nu; p.file.smod=10; p.sw.para=2; p.sw.foldcheck=1; 
p.nc.ilam=1; p.sol.ds=-0.1; p.nc.dsmax=0.1; p.nc.lammin=2; p.sw.bifcheck=2; 
lam=par(1); u=lam*ones(p.np,1); v=(1/lam)*ones(p.np,1); p.u=[u;v;par']; 
p=box2per(p,psw); % p=setfemops(p); % switch on pBC and 
[po,tr,ed]=getpte(p); p.mesh.bp=po; p.mesh.be=ed; p.mesh.bt=tr; p.nc.ngen=1; 
p.plot.bpcmp=7; p.plot.axis='image'; p.plot.cm='cool'; p.usrlam=2:0.1:3.3; 