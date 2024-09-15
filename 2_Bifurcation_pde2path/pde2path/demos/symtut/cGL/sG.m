function r=sG(p,u) % compute pde-part of residual
par=u(p.nu+1:end); f=nodalf(p,u); 
r=par(9)^2*p.mat.K*u(1:p.nu)-p.mat.M0*f-par(9)*par(2)*p.mat.Kx*u(1:p.nu); 