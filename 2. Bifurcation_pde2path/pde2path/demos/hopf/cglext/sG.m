function r=sG(p,u) % compute pde-part of residual
par=u(p.nu+1:end); f=nodalf(p,u); 
r=p.mat.K*u(1:p.nu)-p.mat.M*f-par(6)*p.mat.Kx*u(1:p.nu); 