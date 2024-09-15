function p=Pinit(p,lx,ly,lz, ac, par) % init 
p=stanparam(p); p=stanparamX(p); p.sw.Xfill=1;  p.sw.sfem=-1; 
p.nc.neq=1; p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
% Approximation of Schwatz P surface coming from Fourier series 
% Gandy, et al. "Nodal surface approximations to the PY GY D and I-WP triply periodic 
% minimal surfaces" higher order approximation Schwarz P surface" 
F=@(h,k,l,X,Y,Z) cos(h*X).*(cos(k*Y).*cos(l*Z)+cos(l*Y).*cos(k.*Z))...
    +cos(h*Y).*(cos(k*Z).*cos(l*X)+cos(l*Z).*cos(k*X))...
    +cos(h*Z).*(cos(k*X).*cos(l*Y)+cos(l*X).*cos(k*Y));
p.g=@(p) F(1,0,0,p(:,1),p(:,2),p(:,3))...
    -0.00956335*F(1,1,1,p(:,1),p(:,2),p(:,3))...
    +0.00959441*F(2,1,0,p(:,1),p(:,2),p(:,3))...
    -0.0181364*F(3,0,0,p(:,1),p(:,2),p(:,3))...
    -0.0292538*F(2,2,1,p(:,1),p(:,2),p(:,3));
nit=10; [p.X,p.tri]=dmsurf(p.g,@huniform,ac,[0,0,0;lx,ly,lz],nit);  % nit=2000
% Reflection to full Schwarz P surface
L=[-p.X(:,1) p.X(:,2) p.X(:,3)]; p.X=[p.X;L]; p.tri=[p.tri;p.tri+size(L,1)];
L=[p.X(:,1) -p.X(:,2) p.X(:,3)]; p.X=[p.X;L]; p.tri=[p.tri;p.tri+size(L,1)];
L=[p.X(:,1) p.X(:,2) -p.X(:,3)]; p.X=[p.X;L]; p.tri=[p.tri;p.tri+size(L,1)];
[p.np,~]=size(p.X); p.nu=p.nc.neq*p.np;
[p.X,p.tri]=clean_mesh(p.X,p.tri,'SelfIntersections','ignore'); 
[p.tri,~]=orient_outward(p.X,p.tri,[]);
[p.np,~]=size(p.X); p.nu=p.nc.neq*p.np; [p.nt,~]=size(p.tri); 
u=zeros(p.nu,1); p.u=[u;par];  p.up=p.u; p.sw.orgper=1; p=oosetfemops(p); 
% 3 directional bd-ids 
p.idx1=find(abs(p.X(:,1)-pi)<1e-4 | abs(p.X(:,1)+pi)<1e-4); 
p.idx2=find(abs(p.X(:,2)-pi)<1e-4 | abs(p.X(:,2)+pi)<1e-4);
p.idx3=find(abs(p.X(:,3)-pi)<1e-4 | abs(p.X(:,3)+pi)<1e-4); 
p.idx=[p.idx1:p.idx2;p.idx3]; pde=stanpdeo2D(1,1,5,5); % dummy 
q=p; pde.grid.p=p.X'; pde.grid.t=p.tri'; q.pdeo=pde; 
bcper=[1 2 3]; p=box2per(q,bcper); % periodic BC's in all 3 coordinates 
[Fx,Fy,Fz]=Xfillmat(p); p.Xfillx=Fx; p.Xfilly=Fy; p.Xfillz=Fz; 
m=min(p.X(:,3)); M=max(p.X(:,3)); p.u(p.nu+4)=M-m;  % init parameters
p.u(p.nu+3)=getA(p,p.u); p.u(p.nu+2)=getV(p,p.u);
p.sol.xi=1/p.nu; p.nc.dsmin=1e-4; p.nc.dsmax=1; 
p.nc.dlammax=0.5; p.nc.neig=min(10,p.np); p.sw.para=2; p.sw.jac=0; 
p.sw.spcalc=1; p.plot.pstyle=-1; p.plot.pcmp=1; p.plot.bpcmp=9;  
p.sw.verb=2; p.nc.mu1=0.1; p.nc.mu2=0.01; p.file.smod=1; 
p.sw.bifcheck=2; p.sw.foldcheck=0; p.nc.foldtol=0.1; p.sw.Xfill=1; 
p.plot.auxdict={'H','V','A','\delta','sx','sy','sz'}; 
p.nc.ilam=[4 5 6 7]; p.nc.nq=3; p.fuha.qf=@qf; p.fuha.qfder=@qjac; 
p.sw.jac=0; p.sol.ds=-0.1; p.nc.dsmax=0.1; p.nc.dsmin=1e-6;

