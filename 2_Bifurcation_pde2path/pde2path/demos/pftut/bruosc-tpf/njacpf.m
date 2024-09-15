function [f1u,f1v,f2u,f2v]=njacpf(p,u) % Jacobian 
par=u(p.nu+1:end); u1=u(1:p.np); u2=u(p.np+1:2*p.np);
b=par(2); 
f1u=-(b+1)+2*u1.*u2; 
f1v=u1.^2; 
f2u=b-2*u1.*u2; 
f2v=-u1.^2; 