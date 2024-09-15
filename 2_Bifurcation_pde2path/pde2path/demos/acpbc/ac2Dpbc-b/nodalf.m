function f=nodalf(p,u) % nonlinearity, nodal version
par=u(p.nu+1:end); u1=p.mat.fill*u(1:p.nu); f=par(1)*u1-u1.^3;
