function r=sG(p,u) 
par=u(p.nu+1:end); D=par(5); s=par(6); n=p.np; nu=p.nu; 
K=p.mat.K; Kx=p.mat.Kx; 
u1=u(1:n); u2=u(n+1:2*n); 
bK=[[K 0*K];[0*K D*K]]; 
r=(bK-s*Kx)*u(1:nu)-p.mat.M(1:nu,1:nu)*nodalf(p,u); 