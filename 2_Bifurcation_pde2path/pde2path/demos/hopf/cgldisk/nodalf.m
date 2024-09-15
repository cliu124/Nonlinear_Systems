function f=nodalf(p,u) % nonlinearity for cGL 
par=u(p.nu+1:end); u=p.mat.fill*u(1:p.nu); % extract param., extend sol
r=par(1); nu=par(2); mu=par(3); c3=par(4); c5=par(5); del=par(7);
%nu, pause
n=p.np; u1=u(1:n); u2=u(n+1:2*n); ua=u1.^2+u2.^2; % aux variable |u|^2 
f1=r*u1-nu*u2-ua.*(c3*u1-mu*u2)-c5*ua.^2.*u1; 
f2=r*u2+del^2*nu*u1-ua.*(c3*u2+mu*u1)-c5*ua.^2.*u2; 
f=[f1;f2]; 