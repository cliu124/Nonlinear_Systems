function p=schnakinit(p,n)
p=stanparam(p);screenlayout(p);p.nc.neq=2;p.plot.axis='image';
p.fuha.sG=@sG;p.fuha.sGjac=@sGjac; p.fuha.outfu=@stanbra;
kc=sqrt(sqrt(2)-1);
if n==1;
lx=2*pi/kc;pde=stanpdeo1D(lx,2*lx/100);p.plot.pstyle=1;
end
if n==2
lx=2*pi/kc;ly=2*pi/(sqrt(3)*kc);pde=stanpdeo2D(lx,ly,2*lx/20);
end
if n==3
lx=sqrt(2)*pi/kc/2;pde=stanpdeo3D(lx,lx,lx,2*lx/10);
end
p.Om=(2*lx);
p.np=pde.grid.nPoints;p.nu=p.np*p.nc.neq;p.pdeo=pde;p.sw.sfem=-1;
p=setfemops(p);p.nc.ilam=1;
p.u=[3.3*ones(p.np,1);1/3.3*ones(p.np,1);3.3];
p=setfn(p,['h',num2str(n)]);p.sol.ds=-0.1;p.sol.xi=1/p.nu;p.nc.dsmax=0.1;
p.nc.dsmin=0;p.nc.lammin=0;p.sw.bifcheck=0;p.nc.nsteps=10000;p.file.smod=10;
p.sw.bifcheck=0;p.sw.spcalc=0;
end