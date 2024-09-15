function r=sG(p,u) 
Ks=p.mat.K; M=p.mat.Ms; Kx=p.mat.Kx; 
par=u(p.nu+1:end); al=par(1); eps=par(3); s=par(4); 
n=p.nu/2; u1=u(1:n); 
K=[[-Ks-s*Kx -al*Ks];[-Ks -M]]; F=[-0.5*Kx*(u1.^2)+eps; zeros(n,1)]; 
r=K*u(1:p.nu)-F; 