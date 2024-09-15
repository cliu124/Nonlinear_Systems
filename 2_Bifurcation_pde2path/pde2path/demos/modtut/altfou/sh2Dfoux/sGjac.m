function Gu=sGjac(p,u) % Jacobian 
par=u(p.nu+1:end); u=u(1:p.nu); lam=par(1); c2=par(2); c3=par(3); 
fu=lam+2*c2*u+3*c3*u.^2; 
Fu=spdiags(fu,0,p.nu,p.nu);     
Gu=p.mat.L-Fu; 