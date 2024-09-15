function r=sG(p,u) % compute pde-part of residual
par=u(p.nu+1:end); del=par(7); f=nodalf(p,u); 
%plotsolu(p,u,1,1,1); pause
r=del^2*p.mat.K*u(1:p.nu)-p.mat.M0*f-par(6)*p.mat.Kx*u(1:p.nu); 