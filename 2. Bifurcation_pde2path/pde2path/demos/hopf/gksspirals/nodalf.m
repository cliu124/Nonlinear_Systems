function f=nodalf(p,u) % for rot 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); r=u(p.nu+1); al=u(p.nu+2); 
ua=u1.^2+u2.^2; 
f1=(r+0.5)*u1+u2-ua.*(u1-al*u2); 
f2=r*u2-u1-ua.*(u2+al*u1); 
f=[f1;f2]; 
