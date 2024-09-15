function Gu=sGjac(p,u) % AC, with DBCs
par=u(p.nu+1:end); u=u(1:p.nu); c=par(1); lam=par(2); c2=par(3); c3=par(4); 
fu=lam+2*c2*u+3*c3*u.^2; fu(p.bdi)=0; % local Jac, zeroed on boundary 
n=p.nu; Fu=spdiags(fu,0,n,n); % local Jac, converted to matrix 
Gu=-c*p.mat.L-Fu;     % build Jac 