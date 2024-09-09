function f=nodalf2(p,u) % 'nonlinearity' for cGL (everything but diffusion)  
n=p.nu/2; u1=u(1:n); u2=u(n+1:2*n); par=u(p.nu+1:end); % extract fields 
r=par(1); nu=par(2); mu=par(3); c3=par(4); % and parameters 
ua=u1.^2+u2.^2; % aux variable |u|^2 
fo=fofu2(p,u); 
f1=r*u1-nu*u2-ua.*(c3*u1-mu*u2)+fo; 
f2=r*u2+nu*u1-ua.*(c3*u2+mu*u1)+fo; 
f=[f1;f2]; 