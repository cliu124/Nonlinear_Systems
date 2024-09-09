function p=oosetfemops(p) % Laplacian and D_\phi with NBCs 
[L,Dphi,np,xx,yy,r]=CFLapNBC(p); lx=p.lx; p.r=r(1:p.nr); p.xx=lx*xx; p.yy=lx*yy; 
p.np=np; p.mat.L=sparse(L/(lx^2)); p.mat.M=speye(np); p.mat.Dphi=sparse(Dphi); 