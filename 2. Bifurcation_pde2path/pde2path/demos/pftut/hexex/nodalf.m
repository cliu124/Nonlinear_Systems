function f=nodalf(p,u) % nonlinearity, nodal version
par=u(p.nu+1:end); u=u(1:p.nu); f=par(1)*(u+u.^3); 
