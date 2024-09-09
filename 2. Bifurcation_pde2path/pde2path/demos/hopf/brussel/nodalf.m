function f=nodalf(p,u) % for bru
u1=u(1:p.np); u2=u(p.np+1:2*p.np); u3=u(2*p.np+1:3*p.np); 
par=u(p.nu+1:end); a=par(1); b=par(2); c=par(3); d=par(4); 
f1=a-(1+b)*u1+u1.^2.*u2-c*u1+d*u3; 
f2=b*u1-u1.^2.*u2; 
f3=c*u1-d*u3; 
f=[f1;f2;f3]; 
end