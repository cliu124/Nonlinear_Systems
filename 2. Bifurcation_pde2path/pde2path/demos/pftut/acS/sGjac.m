function Gu=sGjac(p,u)  % generic PDE Jacobian, sfem=+- 1 setting 
par=u(p.nu+1:end); u=u(1:p.nu); R=par(1); lam=par(2); ga=par(3); s=par(4);  
fu=lam+2*ga*u-3*u.^2; 
Fu=spdiags(fu,0,p.nu,p.nu);  
Gu=p.mat.K/R^2-p.mat.M*Fu+s*p.mat.Dphi; 