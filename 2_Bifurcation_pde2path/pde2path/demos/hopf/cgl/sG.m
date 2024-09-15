function r=sG(p,u) % pde-part of residual, generic form, details in nodalf
f=nodalf(p,u); r=p.mat.K*u(1:p.nu)-p.mat.M*f; 