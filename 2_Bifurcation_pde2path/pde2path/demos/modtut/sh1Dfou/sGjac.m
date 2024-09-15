function Gu=sGjac(p,u) % sh1D, (full) Jacobian via F  
n=p.nu; par=u(n+1:end); uf=u(1:n); lam=par(1); c2=par(2); c3=par(3); 
F=p.mat.F; u=F'*uf; fu=lam+2*c2*u+3*c3*u.^2; 
Fu=F*(spdiags(fu,0,n,n)*F');   Gu=p.mat.L-Fu; 