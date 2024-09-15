function r=sG(p,u) % rhs for SH with global coupling in fnl
f=nodalf(p,u)+fnl(p,u); r=p.mat.K*u(1:p.nu)-p.mat.M*f; 