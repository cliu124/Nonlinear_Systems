function p=vgpinit(p) 
% init for vector GP 
p=stanparam(p); screenlayout(p); 
p.nc.neq=4; p.fuha.G=@vgpG; p.fuha.Gjac=@vgpjac; 
p.fuha.outfu=@vgpbra;  p.fuha.postmmod=@gppmm; 
p.file.dircheck=0; pre=sprintf('%s',inputname(1)); p=setfn(p,pre);
lx=5; ly=1*lx; p.nc.tol=1e-6; p.mesh.geo=rec(lx,ly); p.plot.pstyle=3; 
sc=100; qd=[sc, sc, sc, sc]; qm=diag(qd,0); % stiff spring Dirichlet 
bc=gnbc(4,4,qm,zeros(p.nc.neq,1)); 
p.fuha.bc=@(p,u,lam) bc; p.fuha.bcjac=@(p,u) bc;
nx=30;ny=nx;p=stanmesh(p,nx,ny);p=setbmesh(p); p.sol.xi=1/(2*p.nu); 
p.sw.bifcheck=0; p.sw.spcalc=0; p.nc.dsmax=0.1; p.nc.nsteps=30; 
p.plot.pcmp=10; p.nc.imax=50; p.plot.bpcmp=3; p.sol.ds=0.1; 
p.nc.lammax=1.5; p.file.smod=5; p.sw.para=0; p.file.smod=5; 
p.mat.pot=pot(p); p.mat.poti=pdeintrp(p.mesh.p,p.mesh.t,p.mat.pot); % potential
lam=0; pa=1; mu1=2; mu2=2.2; par=[lam;pa;mu1;mu2];  p.nc.ilam=1; % parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; p.r=x.^2+y.^2; p.phi=angle(x+1i*y); 
m=1; n=0; % first dipole (i.e. dipole in u1,v1)
am=pi*factorial(m+1)/2; gm=am; bm=pi*factorial(2*m)/2^(2*m+5); 
dm=pi*factorial(m)/2; a2=(sqrt(dm^2*mu1^2+12*am*gm)+dm*mu1)/(6*gm); 
a=sqrt(a2);A2=(2*a2*gm-dm*mu1)/(bm*3); A=sqrt(A2); 
u=A*tfu(p.r/a,m,n).*cos(m*p.phi); v=0.0*u; u0=[u v]; 
m=1; n=0; % second dipole (i.e. dipole in u2,v2)
am=pi*factorial(m+1)/2; gm=am; bm=pi*factorial(2*m)/2^(2*m+5); 
dm=pi*factorial(m)/2; a2=(sqrt(dm^2*mu2^2+12*am*gm)+dm*mu2)/(6*gm); 
a=sqrt(a2); A2=(2*a2*gm-dm*mu2)/(bm*3); A=sqrt(A2); 
u=A*tfu(p.r/a,m,n).*cos(m*p.phi); v=0.0*u; 
u0=[u0 u v]; p.u=[reshape(u0,p.nc.neq*p.np,1);par]; % append to first and reshape 
plotsol(p,1,1,3); p.tau=0*p.u; 
p.mesh.maxt=2*p.mesh.nt; p.sol.xi=1/p.np; 
res=norm(resi(p,p.u), p.sw.norm); fprintf('initial res=%g\n',res);
p.imax=10;[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res);
plotsol(p,1,1,3);p=meshref(p,'maxt',6000); p=setbmesh(p);

 
