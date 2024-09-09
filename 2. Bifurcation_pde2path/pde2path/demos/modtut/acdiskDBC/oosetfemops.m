function p=oosetfemops(p) % Laplacian and D_\phi with NBCs 
[L,Dphi,np,r,th]=CFLapDBC(p); lx=p.lx; p.r=r(1:p.nr); p.th=th; 
[rr,tt]=meshgrid(p.lx*[p.r],th); [xx,yy]=pol2cart(tt,rr); 
r2=xx.^2+yy.^2; bdi=find(abs(r2-p.lx^2)<1e-6); 
p.bdi=bdi; p.bui=setdiff(1:np+p.na,bdi); 
p.np=np; p.mat.L=sparse(L/(lx^2)); 
p.mat.M=speye(np); p.mat.Dphi=sparse(Dphi); 