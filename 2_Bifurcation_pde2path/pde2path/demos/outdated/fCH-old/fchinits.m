function p=fchinits(p,lx,ly,nx,ny,wi,a,sw,varargin) 
% lx,ly,nx,ny=dom and discret.
% wi(dth), a(mplitude), sw(itch) to describe initial guess
% some other stuff like p.eta1=1; p.eta2=2; p.eps=0.1; p.ga=0; 
% and potential p.wfu=@wfu; set BEFORE init
p=stanparam(p);p.file.dircheck=0;dir=sprintf('%s',inputname(1)); p=setfn(p,dir);
screenlayout(p);
p.nc.neq=2; p.nc.nq=1; % number of add. eq. 
p.fuha.G=@fchG; p.fuha.Gjac=@fchGjac; p.fuha.outfu=@fchbra;% p.fuha.ufu=@fchufu; 
p.nq=1; p.fuha.qf=@fchqf; p.sw.qjac=1; p.fuha.qfder=@fchqjac; 
p.eqn.c=[1;0;0;1;1;0;0;1]; p.eqn.b=0; p.eqn.a=0; p.sw.sfem=1; 
p.fuha.sG=@fchsG;p.fuha.sGjac=@fchsjac;p.fuha.postmmod=@fchpostmeshmod;
p.nc.ilam=1; % primary parameter index 
p.plot.pstyle=2;p.sw.jac=1; p.sw.spcalc=0; p.plot.cm=hot; p.file.smod=5; 
p.sol.ineg=-1; p.plot.bpcmp=3; p.nc.imax=10; p.nc.ngen=1; p.sw.bprint=[1 3]; 
[p.mesh.geo,bc]=recnbc2(lx,ly); p.fuha.bc=@(p,u,lam) bc; p.fuha.bcjac=@(p,u) bc;
p=stanmesh(p,nx,ny); p=setbmesh(p); p.nu=p.nc.neq*p.np; 
p.sol.xi=1/(p.nu); x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)';
p.nc.lammin=-2;p.nc.lammax=4; p.nc.dlammax=1; p.sol.ds=0.1; p.nc.dsmax=0.2;
u=(-1)*ones(p.np,1); v=0*ones(p.np,1); 
if(sw==1) % single stripe 
  u=u+a*1./cosh((x/wi).^1); % set v=eps^2 Del u-W'(u)
  [w,wp,wpp]=p.fuha.wfu(u,p); a=0*p.eps^2;
  v=v+a/wi^2*(2*sinh(x/wi).^2-cosh(x/wi).^2)./(cosh(x/wi).^3)-wp;
end
if(sw==2) % pearl-stripe
   u=u+a*(1+0.2*sin(6*y)).*1./cosh(x/wi).^1; % set v=eps^2 Del u-W'(u)
   %u=u+a*(1+0.2*cos(4*y)).*1./cosh(x/wi).^1;
  [w,wp,wpp]=p.fuha.wfu(u,p); a=0*p.eps^2;
  v=v+a/wi^2*(2*sinh(x/wi).^2-cosh(x/wi).^2)./(cosh(x/wi).^3)-wp;
end
if(sw==3) % meander-stripe
   u=u+a*1./cosh((x-0.3*sin(y/2))/wi).^1; % set v=eps^2 Del u-W'(u)
  [w,wp,wpp]=p.fuha.wfu(u,p); a=0*p.eps^2;
  v=v+a/wi^2*(2*sinh(x/wi).^2-cosh(x/wi).^2)./(cosh(x/wi).^3)-wp;
end
if(sw==5) % quarter circle 
   rp=1.5*lx;r=sqrt((x+lx).^2+(y+ly).^2); 
   u=u+a*1./cosh((r-rp)/wi); 
   [w,wp,wpp]=p.fuha.wfu(u,p); a=p.eps^2; 
   v=v+a/wi^2*(2*sinh((r-rp)/wi).^2-cosh((r-rp)/wi).^2)./cosh((r-rp)/wi).^3; 
   v=v-(r>1)./(r+1).*a/wi.*sinh((r-rp)/wi)./cosh((r-rp)/wi).^2-wp; 
end
if(sw==6) % like quarter square 
   rp=1.4*lx;ep=8; r=((x+lx).^ep+(y+ly).^ep).^(1/ep); 
   u=u+a*1./cosh((r-rp)/wi); 
   [w,wp,wpp]=p.fuha.wfu(u,p); a=p.eps^2; 
   v=v+a/wi^2*(2*sinh((r-rp)/wi).^2-cosh((r-rp)/wi).^2)./cosh((r-rp)/wi).^3; 
   v=v-(r>1)./(r+1).*a/wi.*sinh((r-rp)/wi)./cosh((r-rp)/wi).^2-wp; 
end
if(sw==7) % single spot 
   r=sqrt(x.^2+y.^2); 
   u=u+a*1./cosh(r/wi); 
   [w,wp,wpp]=p.fuha.wfu(u,p); a=p.eps^2; 
   v=0*u; %v+a/wi^2*(2*sinh(r/wi).^2-cosh(r/wi).^2)./cosh(r/wi).^3; 
end
p.u=[u;v; zeros(p.nq,1)]; plotsol(p,1,1,2); pause(1);
p.u(p.nu+1)=p.eta1; % the primary cont.param 
p.u(p.nu+2)=0; % the Lagrange multiplier (gamma)
p.u(p.nu+3)=p.eps; % to later switch to cont. in eps!
p.u(p.nu+5)=p.eta2; 
p=setfemops(p); 
p.nc.ilam=[1 2]; % vector of active aux variables (length p.nq+1) 
p.tau=zeros(length(p.u),1);p.tau(p.nu+1)=1;
p.vol=triint(ones(p.np,1),p.mesh.p,p.mesh.t); 
C=n2triamat(p.mesh.p,p.mesh.t); ta=triar(p.mesh.p,p.mesh.t); % triangle areas 
eta=ta*C; p.eta=[eta zeros(1,p.np)]; mi=p.eta*p.u(1:p.nu)/p.vol; 
fprintf('initial mass=%g, %g\n\n',mi); 
if p.mi~=0; p.u(p.nu+4)=p.mi; else p.u(p.nu+4)=mi; end 
r=resi(p,p.u);res=norm(r,p.sw.norm); fprintf('ini-res=%g\n',res); 
p.sol.xiq=0.05;

