function p=lvinit(p,lx,nx,sw) % init-routine 
p=stanparam(p); screenlayout(p); 
p.nc.neq=4; p.fuha.sG=@lvsG; p.fuha.outfu=@lvbra; % rhs and branch-output
p.fuha.jcf=@lvjcf; p.fuha.con=@lvcon; % current-val-obj and control  
p.usrlam=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8]; 
p.vol=lx; p.xplot=lx; p.sw.spcalc=0; p.sw.jac=0; 
p.file.smod=100; p.plot.bpcmp=1; p.nc.bisecmax=20;
p.nc.dsmin=1e-6; p.nc.dsmax=50; p.nc.lammax=100; p.nc.lammin=-100; p.plot.pstyle=10; 
pde=stanpdeo1D(lx,lx/nx); p.np=pde.grid.nPoints; p.nu=4*p.np; p.pdeo=pde; 
p.sol.xi=1/p.np; po=getpte(p); x=po(1,:)'; lx=max(x); 
p1=20; p2=10; c1=0.1; c2=0.1; d1=1; d2=10; al1=0.4; al2=0.4; beta=0.6; 
ga1=0.1; ga2=0.1; b2=0; b3=0; % b2,b3 unused
h1=1; h2=1-beta; rho=0.03; 
switch sw  % choose initial guess arcoording to switch  
case 1; dir='b1'; sl=0.1/lx; u=h1*ones(p.np,1)+sl*x; v=0.2*h2*ones(p.np,1); 
    lfak=1; l1=50+1*lfak*(ones(p.np,1)-sl*x); l2=10+0.5*lfak*(ones(p.np,1)-sl*x);
    l1=100*(1-x/lx); l2=50*(1-x/lx); % also works, but small imag    
case 2; dir='b2'; sl=1/lx; u=1.25*h1*ones(p.np,1)+sl*x; v=0.1*h2*ones(p.np,1); 
    lfak=2; l1=10+1*lfak*(ones(p.np,1)-sl*x); l2=-10+1*lfak*(ones(p.np,1)-sl*x); % works
    sl=10/lx; lfak=2; l1=10+1*lfak*(ones(p.np,1)-sl*x); l2=-10+1*lfak*(ones(p.np,1)-sl*x); 
end
par=[beta b2 b3 d1 d2 ga1 ga2 rho al1 al2 p1 p2 c1 c2]; p.nc.ilam=1; p.sol.ds=0.05;
%    1     2  3  4  5  6   7   8   9  10  11 12 13 14
p.u=[u;v;l1;l2; par']; 
plotsol(p,1,1,p.plot.pstyle); %return
p.sw.sfem=-1; p=setfemops(p);  % semilin. setting 
r=resi(p,p.u); 
res=norm(r, p.sw.norm); fprintf('initial res=%g\n',res);
p.nc.imax=40; p=setfn(p,dir); %plotsol(p,1,1,10); pause
[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res); plotsol(p,1,1,p.plot.pstyle);
p.nc.imax=10; p.u=real(p.u); 