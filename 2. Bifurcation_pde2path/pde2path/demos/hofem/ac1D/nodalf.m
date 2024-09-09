function f=nodalf(p,u)  % nodal nonlinearity f, called in sG and e2rs
par=u(p.nu+1:end); u=u(1:p.nu); f=par(2)*u+u.^3-par(3)*u.^5; 
