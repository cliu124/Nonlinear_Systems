function f=nodalf(p,u) % nonlinearity, nodal version
par=u(p.nu+1:end); lam=par(2); ga=par(3); 
u=u(1:p.nu); f=lam*u+u.^3-ga*u.^5; 