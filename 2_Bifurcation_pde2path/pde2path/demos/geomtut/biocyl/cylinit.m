function p=cylinit(p,lx,nx,ny,par,sym) % init 
p=stanparam(p); sw.sym=sym; p=stanparamX(p); p.sw.Xcont=2; 
pde=stanpdeo2D(lx,pi,nx,ny,sw); p.sw.jac=0;  % numJac
p.np=pde.grid.nPoints; p.nt=pde.grid.nElements; 
p.fuha.qf=@qfrot; p.fuha.qfder=@qrotder; p.sw.qjac=1;  % for (later) PC 
p.sol.xi=1/p.np; p.nc.neq=2;  p.file.smod=1; p.sw.verb=2; 
p.nc.neig=12; p.nc.eigref=-20; p.sw.bifcheck=2; p.sw.foldcheck=1; 
po=pde.grid.p; x=po(1,:)'; phi=po(2,:)'; %set initial solution
al=par(1); l1=1/(4*al^2); par(2)=l1; uu=0*x+al;
X2=uu.*cos(phi); X3=uu.*sin(phi); p.X=[x, X2, X3]; p.tri=pde.grid.t(1:3,:)'; 
[p.X,p.tri]=clean_mesh(p.X,p.tri,'SelfIntersections','ignore'); 
[p.np,~]=size(p.X); p.nu=2*p.np; [p.nt,~]=size(p.tri); u=zeros(p.np,1);
[q,~,~]=boundary_faces(p.tri); p.idx=unique([q(:,1);q(:,2)]); 
u2=u+0.5/al; % initial H  (guess only) 
p=oosetfemops(p); p.u=[u;u2;par]; 
p.plot.auxdict={'al','l1','c0','s'}; p.plot.bpcmp=5; 
p.fuha.outfu=@hcylbra; p.fuha.e2rs=@e2rsA; p.sw.rlong=1;  
p.nc.mu1=5; p.nc.mu2=0.1; p.nc.bisecmax=5; p.sw.para=1; p.DIR=[]; 
p.nc.foldtol=0.01; p.sw.foldcheck=0; p.nc.ilam=3; p.sol.ds=0.1;