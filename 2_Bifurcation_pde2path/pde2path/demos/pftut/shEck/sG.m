function r=sG(p,u) % rhs for SH, see nodalf
f=nodalf(p,u); l=u(p.nu+3); Ks=p.mat.K; 
K=[0*Ks -l^2*Ks; l^2*Ks p.mat.Ms]; 
r=K*u(1:p.nu)-p.mat.M*f; 