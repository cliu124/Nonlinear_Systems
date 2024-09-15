function r=sG(p,u) % compute pde-part of residual
par=u(p.nu+1:end); f=nodalf(p,u); 
%plotsolu(p,u,1,1,1); pause
r=p.mat.K*u(1:p.nu)-p.mat.M*f-par(6)*p.mat.Krot*u(1:p.nu); 