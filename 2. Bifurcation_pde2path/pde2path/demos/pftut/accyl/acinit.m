function p=acinit(p,lx,ly,nx,par) % AC on cylinder (Om1, pdeo, u1) with 
% lid (Om2, p2, u2). I.e., 2 domains, need to set up 2 PDE objects. Problems 
% coupled via common boundary, which we identify here, and set up useful 
% 'coupling matrices' 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; p.plot.cm='cool'; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.fuha.qf=@qf; p.fuha.qfder=@qfder; 
p.plot.auxdict={'c','lambda','gamma','R','s'}; p.plot.pstyle=-1; % userplot
pde=stanpdeo2D(lx,ly,nx,round(nx*ly/lx)); % domain for cylinder 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; p.sol.xi=1/(p.nu); 
p=box2per(p,1); % pBC on cylinder 
p.nu1=p.nu; p.np1=p.np; % store np and nu for cylinder

[x1,i1]=pde.grid.bdseg(3); t=x1(1,1:end);  % extract (top) boundary, 
% and use as bdry for top disk; i1=indizes of bdry mesh points in 1st domain 
p2.grid=grid2DHU; p2.fem=lagrange12D; 
p2.grid.freeGeometry([cos(t);sin(t)]); 
p.p2=p2; p.np2=p2.grid.nPoints; p.nu2=p.np2; p.nu=p.nu1+p.np2; 
[x2,i2]=p2.grid.bdseg(1);  % i2=indizes of bdry mesh-points of 2nd domain  

p.i1=i1; p.i2=i2; ni=length(i1); % now set up coupling matrices 
Q2=sparse(p.np2,ni); S1=sparse(ni,p.np1); S2=sparse(ni,p.np2); 
for i=1:ni;  Q2(i2(i),i)=1; S1(i,i1(i))=1; S2(i,i2(i))=1; end 
p.Q2=Q2; p.S1=S1; p.S2=S2; p.sf=1e3; 

u=zeros(p.nu1+p.np2,1); p.u=[u; par']; p.nc.ilam=2; p.usrlam=0:2:6; 

p=oosetfemops(p); % now set up system matrices! 
p.nc.nsteps=100; p.sw.foldcheck=0; p.sol.ds=0.1; p.nc.dsmax=0.2; p.nc.lammax=6; 
return
% stuff to test impl. 
x1=getpte(p)'; p.u(1:p.nu1)=cos(x1(:,1)); 
x2=p2.grid.p; x2=x2'; phi=angle(x2(:,1)+1i*x2(:,2)); 
p.u(p.nu1+1:p.nu1+p.np2)=cos(phi); plotsol(p,1,1,10); 



