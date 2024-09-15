function Gu=sGjac(p,u) % Jacobian 
n=p.nu; par=u(n+1:end); uf=u(1:n); c=par(1); lam=par(2); c2=par(3); c3=par(4); 
F=p.mat.F; u=F'*uf; fu=lam+2*c2*u+3*c3*u.^2; 
Fu=F*(spdiags(fu,0,n,n)*F'); % local Jac, turned to matrix 
Gu=c*spdiags(p.mat.k2,0,n,n)-Fu; 