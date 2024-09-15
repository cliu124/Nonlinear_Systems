function r=sG(p,u) % compute pde-part of residual
par=u(p.nu+1:end); Ks=p.mat.K; 
K=[[par(1)*Ks par(2)*Ks];[par(3)*Ks Ks]]; 
f=nodalf(p,u); 
r=K*u(1:p.nu)-p.mat.M*f; 