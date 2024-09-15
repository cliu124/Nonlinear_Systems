function p=nodinitl(p,lx,ly,nx,ny,par,sym) % init 
p=stanparam(p); sw.sym=sym; pde=stanpdeo2D(lx,ly,nx,ny,sw); 
p.fuha.outfu=@cmcbra; p.fuha.sG=@sGnodD; p.fuha.sGjac=@sGnodDjac; 
p.np=pde.grid.nPoints; p.nu=p.np; p.nt=pde.grid.nElements; %store dimensions
p.sol.xi=1/p.nu; p.nc.neq=1; p.fuha.e2rs=@e2rs; p.sw.sfem=-1;
p.sw.spcalc=1; p.plot.pstyle=-1; p.nc.neig=10; p.file.smod=1; p.sw.verb=2; p.nc.lammax=2.15; 
p.nc.lammin=-10^3; p.sw.bifcheck=2; p.sw.foldcheck=1; p.sf=1e0; p.sw.Xcont=1; 
%set initial solution
po=pde.grid.p; x=po(1,:)'; phi=po(2,:)'; p.tri=pde.grid.t(1:3,:)'; 
R=1; r=0.2; x1=min(x); x2=max(x); 
xs=2; ns=20; del=0.1; % extend x (axial) a bit to adjust the parametrization 
xa=linspace(xs*x1, x1-del,ns); xb=linspace(x2+del,xs*x2,ns); xp=[xa'; (x); xb']; 
k2=(R^2-r^2)/(R^2); k=sqrt(k2); 
xp=(x); ns=0; e1=ellipticE(xp,k); f1=ellipticF(xp,k);
d=sqrt(1-k2*sin(xp).^2); x1=1*R*e1-r*f1-R*k2*sin(xp).*cos(xp)./d; 
r1=r./d;   % nodary following Mladenov 2002 
x2=x1(ns+1:end-ns); r2=r1(ns+1:end-ns); % chop off extensions again: 
par(3)=r2(1); par(1)=1./(R-r); % store H 
p.X=[r2.*cos(phi), r2.*sin(phi),x2];p.X=[r2.*cos(phi), r2.*sin(phi),x2]; 
[p.X,p.tri]=clean_mesh(p.X,p.tri,'SelfIntersections','ignore'); 
p.np=size(p.X,1); p.nu=p.np; p.nt=size(p.tri,1); u=zeros(p.np,1);
[p.DIR,~,~]=boundary_faces(p.tri); p.idx=unique([p.DIR(:,1);p.DIR(:,2)]);
p=oosetfemops(p); p.u=[u;par]; p.u(p.np+2)=getV(p,p.u); p.u(p.np+3)=getA(p,p.u); 
p.nc.del=1e-2; p.plot.auxdict={'H','V','A','r','l'}; 