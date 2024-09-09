function p=nodbuckinit(p,lx,ly,nx,ny,par,sym) % init 
p=stanparam(p); p=stanparamX(p); sw.sym=sym; pde=stanpdeo2D(lx,ly,nx,ny,sw); 
p.fuha.sG=@sGnodpBC; p.fuha.sGjac=@sGnodpBCjac;
p.np=pde.grid.nPoints; p.nu=p.np; p.nt=pde.grid.nElements; %store dimensions
p.sol.xi=1/p.nu; p.nc.neq=1; p.fuha.e2rs=@e2rs; p.sw.sfem=-1;
p.sw.spcalc=1; p.plot.pstyle=-1; p.nc.neig=20; p.file.smod=1; p.sw.verb=2; 
p.sw.bifcheck=2; p.sw.foldcheck=0; 
%set initial solution, again using KPP-parametr.
po=pde.grid.p; p.tri=pde.grid.t(1:3,:)'; h=par(1); a=par(4); x=po(1,:)'+lx; 
phi=po(2,:)'; x1=(cos(x)+sqrt(a+cos(x).^2))/(2*abs(h)); 
l=size(x,1); x2=0*x1; % zeros(size(x)); 
for i=1:l % z=z(x) by (numerical) integration 
  y=@(t) 1./(2*abs(h)).*(cos(t)+sqrt(cos(t).^2+a)).*cos(t)./sqrt(cos(t).^2+a);
  x2(i)=integral(y,0,x(i));
end
p.X=[x1.*cos(phi), x1.*sin(phi),x2]; % rotate 
m=min(p.X(:,3)); M=max(p.X(:,3)); par(6)=M-m; % initial height
par(5)=min(sqrt(p.X(:,1).^2+p.X(:,2).^2)); % initial inner radius 
p.X=[p.X(:,1) p.X(:,2) p.X(:,3)-par(6)/2];
[p.X,p.tri]=clean_mesh(p.X,p.tri,'SelfIntersections','ignore'); 
p.np=size(p.X,1); p.nu=p.np; p.nt=size(p.tri,1); u=zeros(p.np,1);
bdf=boundary_faces(p.tri); p.idx=unique([bdf(:,1);bdf(:,2)]);
p=oosetfemops(p); p.u=[u;par]; 
% pdeo etc only needed once, so here make dummy q, later discarded 
q=p; q.pdeo=pde; q.pdeo.grid.p=q.X'; q.pdeo.grid.t=p.tri; 
p=box2per(q,3); % generate (initial) p.mat.fill and p.mat.drop 
p.u(p.nu+2)=getV(p,p.u); p.nc.bisecmax=10; 
p.plot.auxdict={'H','V','A','r','l','delta'};  p.nc.tol=1e-6; 
p.plot.bpcmp=19; p.sw.jac=0; p.nc.dsmax=0.04; p.sw.para=2; 
p.nc.mu1=0.25; p.nc.mu2=0.1; p.nc.foldtol=1e-2; p.nc.bisecmax=10;
