function r=sG(p,u)  % PDE rhs 
n=p.nu; par=u(n+1:end); u=u(1:n); c=par(1); lam=par(2); c2=par(3); c3=par(4); 
f=lam*u+c2*u.^2+c3*u.^3; f(p.bdi)=0; % zero nonlin on bdry 
r=-c*p.mat.L*u-f; 