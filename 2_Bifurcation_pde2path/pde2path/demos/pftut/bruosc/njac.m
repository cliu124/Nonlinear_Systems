function [f1u,f1v,f2u,f2v]=njac(p,u) % Jacobian 
par=u(p.nu+1:end); u1=u(1:p.nu/2); u2=u(p.nu/2+1:p.nu);
a=par(1); b=par(2); 
f1u=-(b+1)+2*u1.*u2; 
f1v=u1.^2; 
f2u=b-2*u1.*u2; 
f2v=-u1.^2; 