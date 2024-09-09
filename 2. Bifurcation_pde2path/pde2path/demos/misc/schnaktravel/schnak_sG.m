function r=schnak_sG(p,u) % nodal form of residual
f=nodalf(p,u); par=u(p.nu+1:end); 
K=par(3)^2*p.mat.K-par(3)*par(2)*p.mat.Kx; 
r=K*u(1:p.nu)-p.mat.M*f; 



