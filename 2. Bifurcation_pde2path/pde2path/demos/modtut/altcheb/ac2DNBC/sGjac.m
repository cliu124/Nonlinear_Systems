function Gu=sGjac(p,u) % AC 
par=u(p.nu+1:end); u=u(1:p.nu); c=par(1); lam=par(2); c2=par(3); c3=par(4); n=p.nu; 
fu=lam+2*c2*u+3*c3*u.^2; Fu=spdiags(fu,0,n,n); % loca Jac, converted to matrix 
Gu=-c*p.mat.L-Fu; % build Jac only for bulk (u=0 fixed at boundary) 