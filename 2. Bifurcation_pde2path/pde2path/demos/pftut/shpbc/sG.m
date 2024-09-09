function r=sG(p,u) % rhs for SH, see nodalf
par=u(p.nu+1:end); sx=par(3); sy=par(4); 
f=nodalf(p,u); u=u(1:p.nu); 
r=p.mat.K*u-p.mat.M*f+sx*p.mat.Dx*u+sy*p.mat.Dy*u; 