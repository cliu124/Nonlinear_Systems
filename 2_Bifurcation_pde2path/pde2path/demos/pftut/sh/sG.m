function r=sG(p,u) % rhs for SH, see nodalf
f=nodalf(p,u); r=p.mat.K*u(1:p.nu)-p.mat.M*f; 