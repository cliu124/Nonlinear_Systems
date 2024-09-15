function r=sG(p,u) % Schnakenberg 
f=nodalf(p,u); par=u(p.nu+1:end); % lam,sig,d 
K=kron([[1,0];[0,par(3)]],p.mat.K); 
r=K*u(1:p.nu)-p.mat.M*f; 