function q=qf(p,u) % 3 translations 
u=u(1:p.nu); N=getN(p,p.X); 
grX=grad(p.X,p.tri); grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:); grXz=grX(2*p.nt+1:3*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; grXz=E*grXz;
dphix=grXx*p.X; dphix=dot(dphix,N,2); dphiy=grXy*p.X; dphiy=dot(dphiy,N,2);
dphiz=grXz*p.X; dphiz=dot(dphiz,N,2);
qx=(p.mat.fill'*dphix)'*u; qy=(p.mat.fill'*dphiy)'*u; qz=(p.mat.fill'*dphiz)'*u;
q=[qx; qy; qz]; 