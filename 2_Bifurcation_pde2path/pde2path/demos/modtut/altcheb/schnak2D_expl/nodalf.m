function f=nodalf(p,u) % for Schnakenberg 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end); % lam,sig,d 
f1=-u1+u1.^2.*u2+par(2)*(u1-u2.^(-1)).^2; 
f2=par(1)-u1.^2.*u2-par(2)*(u1-u2.^(-1)).^2; 
f=[f1; f2]; 