function [f1u,f1v,f2u,f2v]=njac(p,u) % local (no spat.derivatives) Jacobian 
par=u(p.nu+1:end); n=p.np; u=p.mat.fill*u(1:p.nu); 
r=par(1); nu=par(3); mu=par(4); c3=par(5); c5=par(6); % and parameters 
u1=u(1:n); u2=u(n+1:2*n); ua=u1.^2+u2.^2; 
f1u=r-2*u1.*(c3*u1-mu*u2)-c3*ua-4*c5*ua.*u1.^2-c5*ua.^2; 
f1v=-nu-2*u2.*(c3*u1-mu*u2)+mu*ua-4*c5*ua.*u1.*u2; 
f2u=nu-2*u1.*(c3*u2+mu*u1)-mu*ua-4*c5*ua.*u1.*u2; 
f2v=r-2*u2.*(c3*u2+mu*u1)-c3*ua-4*c5*ua.*u2.^2-c5*ua.^2;  
end
