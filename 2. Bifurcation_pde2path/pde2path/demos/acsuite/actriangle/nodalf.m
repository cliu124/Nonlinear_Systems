function f=nodalf(p,u) % nonlinearity, nodal version
par=u(p.nu+1:end); u=u(1:p.nu); f=par(2)*u+u.^3-par(3)*u.^5; 
