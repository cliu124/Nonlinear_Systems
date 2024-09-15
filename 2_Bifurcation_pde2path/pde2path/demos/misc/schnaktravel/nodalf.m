function f=nodalf(p,u) % schnak-nonlin. nodal version
np=p.nu/2; u1=u(1:np); u2=u(np+1:2*np); par=u(p.nu+1:end);  
f1=-u1+u1.^2.*u2; f2=par(1)-u1.^2.*u2; f=[f1;f2];
