function p=gpinit(p) 
% init-routine 
p=stanparam(p); screenlayout(p); 
p.nc.neq=2; p.fuha.G=@gpG; p.fuha.Gjac=@gpjac; 
%p.fuha.sG=@gpsG; p.fuha.sGjac=@gpsjac; 
p.fuha.outfu=@gpbra; p.fuha.postmmod=@gppmm; 
p.file.dircheck=0; pre=sprintf('%s',inputname(1)); p=setfn(p,pre);
lx=5; ly=1*lx; p.mesh.geo=rec(lx,ly); p.plot.cm='cool'; p.plot.pstyle=3; 
bc=gnbc(2,4,[[100 0];[0 100]],[0;0]); p.fuha.bc=@(p,u,lam) bc; p.fuha.bcjac=@(p,u) bc; 
nx=30; ny=nx; p=stanmesh(p,nx,ny); p.sol.xi=1/(2*p.nu); p=setbmesh(p); 
p.sw.bifcheck=0; p.sw.spcalc=0; 
p.plot.pcmp=10; p.nc.imax=50; p.plot.bpcmp=3; 
p.nc.lammax=1.5; p.file.smod=5; p.sw.para=0; p.file.smod=5; 
p.mat.pot=pot(p); p.mat.poti=pdeintrp(p.mesh.p,p.mesh.t,p.mat.pot); % potential
lam=0; pa=1; mu=2; par=[lam;pa;mu]; % parameters 
p.sol.ds=0.1; p.sol.xi=1/p.np;
x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; r=x.^2+y.^2; 
u=2./cosh(r); % monopole (to work with m=1,2 first meshref the monopole) 
v=0.0*u; p.u=[u(:); v(:); par]; p.nc.ilam=1;
p.tau=zeros(p.nu+1,1); p.tau(p.nu+1)=1; 
res=norm(resi(p,p.u), p.sw.norm); fprintf('initial res=%g\n',res); 
p.tau=zeros(length(p.u),1); p.tau(length(p.u))=1;
p.jsw=1; p.imax=10; [p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res);
plotsol(p,1,1,3); p.tau=0*p.u; p=meshref(p,'maxt',2500); 
% now initial quadpole guess on the already refined mesh 
x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; p.r=x.^2+y.^2; p.phi=angle(x+1i*y); 
m=2; n=0; p.lam=0;
am=pi*factorial(m+1)/2; gm=am; bm=pi*factorial(2*m)/2^(2*m+5); 
dm=pi*factorial(m)/2; a2=(sqrt(dm^2*mu^2+12*am*gm)+dm*mu)/(6*gm); 
a=sqrt(a2)*(n+1); A2=(2*a2*gm-dm*mu)/(bm*3); A=sqrt(A2); 
u=A*tfu(p.r/a,m,n).*(cos(m*p.phi))*(n+1); v=0.0*u; u0=[u v]; 
p.u=[u(:);v(:); par];  plotsol(p,1,1,3); p.tau=0*p.u; 
p.nc.ngen=10;p.mesh.maxt=2*p.mesh.nt; p.sol.xi=1/p.np; 
res=norm(resi(p,p.u), p.sw.norm); fprintf('initial res=%g\n',res);
p.imax=10;[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res);
plotsol(p,1,1,3); p=meshref(p,'maxt',5000); p=setbmesh(p);


