function r=sG(p,u)  % PDE rhs 
par=u(p.nu+1:end); u=u(1:p.nu); lam=par(1); c2=par(2); c3=par(3); % split 
f=lam*u+c2*u.^2+c3*u.^3; % "nonlinearity"  
r=p.mat.L*u-f;     % residual 
  