function f=nodalf(p,u)  % nodal nonlinearity f, called in sG and e2rs
par=u(p.nu+1:end); u=u(1:p.nu); lam=par(1); c2=par(2); c3=par(3); 
f=lam*u+c2*u.^2+c3*u.^3; 
