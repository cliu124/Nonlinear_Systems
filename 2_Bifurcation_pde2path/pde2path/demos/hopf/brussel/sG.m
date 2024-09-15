function r=sG(p,u) % compute pde-part of residual
f=nodalf(p,u); r=p.mat.K*u(1:p.nu)-p.mat.M*f; 