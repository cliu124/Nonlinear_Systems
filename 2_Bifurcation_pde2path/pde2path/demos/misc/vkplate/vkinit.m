function p=vkinit(p) 
% von karman init-routine, partially clamped plate 
p=stanparam(p); screenlayout(p); p.nc.neq=10; p.fuha.G=@vkf; p.fuha.Gjac=@vkjac; 
p.fuha.outfu=@stanbra; p.fuha.sG=@vksG; p.fuha.sGjac=@vksGjac; 
ly=pi/2; lx=1.6*ly; p.nc.tol=1e-8; p.mesh.geo=rec(lx,ly); 
sc=10^3; zv=zeros(p.nc.neq,1); % define BC: 
% III(u), I(v), on x=0,1, uyy=0; (vert. boundaries) 
% III(u), I(v), on y=0,1, uxx=0; (horiz. boundaries) 
qdv=[0 0 sc sc 0 sc 0 sc sc 0]; 
qdh=[sc sc sc sc sc sc 0 sc sc 0]; 
% force u=0 on vertical boundaries
qmv=diag(qdv,0); qmh=diag(qdh,0); qmv(2,1)=sc; 
% ordering of boundaries: bottom, right, top, left
bc=gnbc(10,qmh,zv,qmv,zv,qmh,zv,qmv,zv);
p.fuha.bc=@(p,u,lam) bc; p.fuha.bcjac=@(p,u,lam) bc; 
nx=25;ny=round(nx*ly/lx);p=stanmesh(p,nx,ny); p=setbmesh(p); 
p.file.dirchecksw=0; dir=sprintf('%s',inputname(1)); p=setfn(p,dir);
p.nc.dsmin=0.001; p.nc.dsmax=1; p.nc.dlammax=2; p.nc.nsteps=10; p.nc.lammax=100; 
p.nc.amod=0; p.sw.bifcheck=1; p.sw.spcalc=0; % 0/1
p.plot.pcmp=1; p.nc.imax=10; p.plot.pstyle=2; p.file.smod=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%
lam=5; p.sol.ds=0.1; p.sol.xi=1/p.np; del=0.05; 
u=zeros(1,p.np); u0=[u u u u u u u u u u]; 
p.u=[reshape(u0,p.nc.neq*p.np,1); lam; del]; p.nc.ilam=1;  
plotsol(p,1,1,p.plot.pstyle); p.tau=0*p.u; 
%%% sfem-setup 
c=zeros(400,1); 
% first column of C as written in manual 
c(1)=1; c(4)=1; c(5)=0; c(17)=1; c(24)=1; c(26)=0.5;c(27)=0.5; 
c(45)=1; c(48)=1; % second column 
c(89)=1; c(92)=1; c(109)=1; c(116)=1; c(118)=0.5; c(119)=0.5; % third
c(133)=1; c(136)=1; c(177)=del; c(180)=del; % 4th 
c(221)=del; c(224)=del; c(265)=del; c(268)=del; c(309)=del; c(312)=del; 
c(353)=del; c(356)=del; c(397)=del; c(400)=del;   % to 10th 
a=zeros(100,1); a(11)=1; a(33)=1; a(45)=1; a(56)=1; 
a(67)=1; a(78)=1; a(89)=1; a(100)=1;
p.eqn.c=c; p.eqn.a=a; p.eqn.b=0; 
p.eqn.c2=zeros(400,1); p.eqn.c2(5)=1; 
p.sw.sfem=1; p=setfemops(p); %p.sw.sfem=0; 
