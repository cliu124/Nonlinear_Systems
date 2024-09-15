function f=nodalf(p,u) % for Schnakenberg 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end); 
a=par(1); b=par(2); 
f1=a-(b+1)*u1+u1.^2.*u2;
f2=b*u1-u1.^2.*u2; 
f=[f1; f2]; 