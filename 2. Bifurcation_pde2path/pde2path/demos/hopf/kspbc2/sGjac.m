function Gu=sGjac(p,u)
Ks=p.mat.K; M=p.mat.Ms; Kx=p.mat.Kx; par=p.u(p.nu+1:end); al=par(1); s=par(4); 
n=p.nu/2; u1=u(1:n);
K=[[-Ks-s*Kx -al*Ks];[-Ks -M]]; 
%F1u=-0.5*(spdiags(u1,0,n,n)*Kx+spdiags(Kx*u1,0,n,n));
F1u=-Kx*spdiags(u1,0,n,n); 
Fu=[[F1u 0*F1u];[0*F1u 0*F1u]]; 
Gu=K-Fu; 