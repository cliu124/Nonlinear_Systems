function [f1u,f1v,f2u,f2v]=njac(p,u) % Jacobian for Schnakenberg
%u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
par=u(p.nu+1:end); u1=u(1:p.nu/2); u2=u(p.nu/2+1:p.nu);
f1u=-1+2*u1.*u2+2*par(2)*(u1-u2.^(-1));
f1v=u1.^2+2*par(2)*u2.^(-2).*(u1-u2.^(-1)); 
f2u=-2*u1.*u2-2*par(2)*(u1-u2.^(-1)); 
f2v=-u1.^2-2*par(2)*u2.^(-2).*(u1-u2.^(-1)); 