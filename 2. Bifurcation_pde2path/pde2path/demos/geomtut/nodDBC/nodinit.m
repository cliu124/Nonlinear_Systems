function p=nodinit(p,t,ly,nx,ny,par,sym,varargin) % init, 'short' nodoids 
% t=l for cylinder, flagged by icsw=0 (default for varargin=[]
try icsw=varargin{1}; catch icsw=0; end 
p=stanparam(p); p=stanparamX(p); sw.sym=sym; pde=stanpdeo2D(t,ly,nx,ny,sw); 
p.fuha.outfu=@cmcbra; p.fuha.sG=@sGnodD; p.fuha.sGjac=@sGnodDjac; 
p.nc.ilam=[3 1]; p.nc.nq=1; p.fuha.qf=@qfA; p.fuha.qfder=@qjacA; p.sw.qjac=1; 
p.np=pde.grid.nPoints; p.nu=p.np; p.nt=pde.grid.nElements; %store dimensions
p.sol.xi=1/p.nu; p.fuha.e2rs=@e2rs; p.plot.pstyle=-1; p.plot.bpcmp=1; 
p.nc.neig=20; p.file.smod=1; p.sw.verb=2; p.nc.lammax=100; p.nc.lammin=0; 
p.sw.bifcheck=2; p.sw.foldcheck=1; p.sf=1e0; p.tri=pde.grid.t(1:3,:)';  
% set initial solution
h=par(1); a=par(4); po=pde.grid.p; x=po(1,:)'; phi=po(2,:)'; 
if icsw==0; p.X=[a.*cos(phi), a.*sin(phi),x]; % just cylinder 
else % KPP17 parametrization of nodoids 
  x1=(cos(x)+sqrt(a+cos(x).^2))/(2*abs(h)); x2=zeros(size(x)); xl=size(x,1); 
  for i=1:xl 
    y=@(t) 1./(2*abs(h)).*(cos(t)+sqrt(cos(t).^2+a)).*cos(t)./sqrt(cos(t).^2+a);
    x2(i)=integral(y,0,x(i));  
  end
  p.X=[x1.*cos(phi), x1.*sin(phi), x2]; % initial X 
end
par(7)=max(p.X(:,3))-min(p.X(:,3)); % needed in getV 
[p.X,p.tri]=clean_mesh(p.X,p.tri,'SelfIntersections','ignore'); 
p.np=size(p.X,1); p.nu=p.np; p.nt=size(p.tri,1); u=zeros(p.np,1);
[p.DIR,~,~]=boundary_faces(p.tri); p.idx=unique([p.DIR(:,1);p.DIR(:,2)]);
p=oosetfemops(p); p.u=[u;par]; p.u(p.np+2)=getV(p,p.u); p.u(p.np+3)=getA(p,p.u); 
p.plot.auxdict={'H','V','A','r','l'}; p.nc.foldtol=0.01; 
p.nc.dlammax=5; p.sw.jac=1; p.sol.ds=1; p.nc.dsmax=5; p.nc.mu2=0.01; p.nc.mu1=1; 