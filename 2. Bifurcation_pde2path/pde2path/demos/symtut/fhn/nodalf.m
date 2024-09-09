function f=nodalf(p,u) % nodal version of nonlinearity
par=u(p.nu+1:end); u1=u(1:p.nu/p.nc.neq); u2=u(p.nu/p.nc.neq+1:p.nu); 
f1=u1-u1.^3 - par(1)*(par(3) + par(4)*u2 + par(5)*u2.^2 + par(6)*u2.^3); 
f2=par(1)^2*(u1-u2);
f=[f1;f2];
end