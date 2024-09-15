function Gu=sGjac(p,u)
K=p.mat.K; M0=p.mat.M0; Kx=p.mat.Kx; par=u(p.nu+1:end); al=par(1); s=par(4); 
n=p.nu; u=u(1:n);
Fu=M0*(Kx*spdiags(u,0,n,n)); 
Gu=al*K^2-M0*K+Fu+s*M0*Kx; 