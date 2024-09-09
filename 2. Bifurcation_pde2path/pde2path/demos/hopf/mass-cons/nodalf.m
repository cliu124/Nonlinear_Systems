function f=nodalf(p,u) % nonlinearity for cGL 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end); % extract fields 
a1=par(5); b1=par(6); % and parameters 
f1=a1*u1-u1.^3+b1*u1.*u2;  f2=-f1;  f=[f1;f2]; 