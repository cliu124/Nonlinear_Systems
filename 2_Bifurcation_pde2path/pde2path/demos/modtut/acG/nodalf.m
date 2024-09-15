function f=nodalf(p,u)  % nodal nonlinearity f, called in sG 
par=u(p.nu+1:end); u=u(1:p.nu); lam=par(2); c2=par(3); c3=par(4); 
f=lam*u+c2*u.^2+c3*u.^3; 
