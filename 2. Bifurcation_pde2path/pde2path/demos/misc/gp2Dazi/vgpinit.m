function p=vgpinit(p,lx,nx,par,sw,m) 
p=stanparam(p); screenlayout(p); p.nc.neq=2; 
p.fuha.sG=@gpsG; p.fuha.sGjac=@gpsGjac; p.sw.sfem=-1; p.fuha.outfu=@gpbra; 
p.fuha.qf=@qf; p.fuha.qfder=@qfder; p.sw.qjac=1; 
pde=stanpdeo2D(lx,lx,nx,nx,sw); p.plot.cm='cool'; p.plot.pstyle=-1; 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=2*p.np; p.sol.xi=1/(p.nu); 
p.plot.pcmp=10; p.plot.bpcmp=0; p.nc.lammax=1; p.file.smod=2; 
p.mat.pot=pot(p); % potential 
p.sol.ds=0.1; p.sol.xi=1/p.np; p.sw.bifcheck=0; p.nc.dsmax=0.1; 
p=oosetfemops(p); po=getpte(p); x=po(1,:)'; y=po(2,:)'; r=x.^2+y.^2; 
u=2./cosh(r); % monopole (to work with m=1,2 first meshref the monopole) 
v=0.0*u; p.u=[u(:); v(:); par']; p.nc.ilam=1; 
p.tau=zeros(p.nu+1,1); p.tau(p.nu+1)=1; p.sf=1e3; 
p.nc.ngen=1; p.nc.maxt=20000; p.nref=500; p.fuha.e2rs=@e2rs_rad; 
% now set vector valued and generate iguess 
u0=p.u(1:p.nu); p.nc.neq=4; p.nu=p.nc.neq*p.np; 
po=getpte(p); x=po(1,:)'; y=po(2,:)'; 
p.r=x.^2+y.^2; p.phi=angle(x+1i*y); 
n=0; mu1=par(2); mu2=par(5);  
p.fuha.sG=@vgpsG; p.sw.jac=0; % p.fuha.sGjac=@vgpsGjac; % still to do! 
p.sw.sfem=-1; p.fuha.outfu=@vgpbra; 
p.fuha.qf=@vqf; p.fuha.qfder=@vqfder;
M=p.mat.M; p.mat.M=[M 0*M; 0*M M]; 
am=pi*factorial(m+1)/2; gm=am; bm=pi*factorial(2*m)/2^(2*m+5); 
dm=pi*factorial(m)/2; a2=(sqrt(dm^2*mu1^2+12*am*gm)+dm*mu1)/(6*gm); 
a=sqrt(a2);A2=(2*a2*gm-dm*mu1)/(bm*3); A=sqrt(A2); 
u=A*tfu(p.r/a,m,n).*cos(m*p.phi); v=0.0*u; 
n=1; % second dipole (i.e. dipole in u2,v2)
am=pi*factorial(m+1)/2; gm=am; bm=pi*factorial(2*m)/2^(2*m+5); 
dm=pi*factorial(m)/2; a2=(sqrt(dm^2*mu2^2+12*am*gm)+dm*mu2)/(6*gm); 
a=sqrt(a2); A2=(2*a2*gm-dm*mu2)/(bm*3); A=sqrt(A2); 
u=A*tfu(p.r/a,m,n).*cos(m*p.phi); v=0.0*u; 
p.u=[u0;u;v;par'];  plotsol(p,1,3,3); 
r=resi(p,p.u); plotsol(p,1,1,3); p.tau=0*p.u; 
res=norm(r, p.sw.norm); fprintf('initial res=%g\n',res); 
[p.u,res]=nloop(p,p.u); fprintf('first res=%g\n',res);
plotsol(p,1,1,3); plotsol(p,6,3,3); 


