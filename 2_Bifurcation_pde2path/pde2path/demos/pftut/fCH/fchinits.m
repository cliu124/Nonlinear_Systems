function p=fchinits(p,lx,ly,nx,ny,par,wi,a,sw,varargin) 
% lx,ly,nx,ny=dom and discret.
% wi(dth), a(mplitude), sw(itch) to describe initial guess 
p=stanparam(p); screenlayout(p); p.nc.neq=2; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.outfu=@fchbra; p.fuha.e2rs=@e2rs; 
p.nc.nq=2; p.fuha.qf=@fchqf2; p.sw.qjac=1; p.fuha.qfder=@fchqjac2; 
p.sw.sfem=-1; p.nc.neig=20; p.sw.bifcheck=1; p.nc.nsteps=100; 
p.plot.pstyle=2; p.sw.jac=1; p.sw.spcalc=1; p.plot.cm=hot; p.file.smod=5; 
p.plot.bpcmp=3; p.nc.ngen=1; p.sw.bprint=[1 3]; 
sym.sym=2; pde=stanpdeo2D(lx,ly,nx,ny,sym); p.Om=4*lx*ly; 
p.np=pde.grid.nPoints; p.pdeo=pde; p.nu=p.np*p.nc.neq; 
xp=getpte(p);  x=xp(1,:)'; y=xp(2,:)';
p.nc.lammin=0; p.nc.lammax=4; p.nc.dlammax=1; p.sol.ds=0.1; p.nc.dsmax=0.2;
epsi=par(3); u=(-1)*ones(p.np,1); v=0*ones(p.np,1); 
switch sw % select initial guess 
  case 1; u=u+a*1./cosh((x/wi).^1); % straight channel 
   [w,wp,wpp]=p.fuha.wfu(u,p); a=0*epsi^2;
   v=v+a/wi^2*(2*sinh(x/wi).^2-cosh(x/wi).^2)./(cosh(x/wi).^3)-wp;
   p.plot.axis='tight';
  case 2; rp=0.5*lx; r=sqrt(x.^2+(y+ly).^2); % half circle   
   u=u+a*1./cosh((r-rp)/wi); [w,wp,wpp]=p.fuha.wfu(u,p); a=epsi^2; 
   v=v+a/wi^2*(2*sinh((r-rp)/wi).^2-cosh((r-rp)/wi).^2)./cosh((r-rp)/wi).^3; 
   v=v-(r>1)./(r+1).*a/wi.*sinh((r-rp)/wi)./cosh((r-rp)/wi).^2-wp; 
   p.plot.axis='image';
  case 3;  r=sqrt(x.^2+y.^2); u=u+a*1./cosh(r/wi);  % single spot 
   [w,wp,wpp]=p.fuha.wfu(u,p); a=epsi^2; v=0*u; 
  case 4; rp=1.5*lx;r=sqrt((x+1.25*lx).^2+(y+ly).^2); % distorted quarter circle   
   u=u+a*1./cosh((r-rp)/wi); [w,wp,wpp]=p.fuha.wfu(u,p); a=epsi^2; 
   v=v+a/wi^2*(2*sinh((r-rp)/wi).^2-cosh((r-rp)/wi).^2)./cosh((r-rp)/wi).^3; 
   v=v-(r>1)./(r+1).*a/wi.*sinh((r-rp)/wi)./cosh((r-rp)/wi).^2-wp; 
end
p.u=[u;v; par];  plotsol(p,1,1,2); p=oosetfemops(p); 
p.tau(1:p.nu+2)=0*p.u(1:p.nu+2); p.nc.ilam=[1 2 6]; p.sol.xiq=0.05;
p.vol=4*lx*ly; mi=sum(p.mat.M*p.u(1:p.nu))/p.vol; 
fprintf('initial mass=%g, %g\n\n',mi); ium=par(4); % initial mass given by user 
if ium~=0; p.u(p.nu+4)=ium; else p.u(p.nu+4)=mi; end 
r=resi(p,p.u); res=norm(r,p.sw.norm); fprintf('ini-res=%g\n',res); 