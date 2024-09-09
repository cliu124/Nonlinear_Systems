function Gu=sGjac(p,u) % matrix free Jacobian, 
global p2pglob; % p2pglob.fu further used in lssgmres 
n=p.nu; par=u(n+1:end); uf=u(1:n); lam=par(1); c2=par(2); c3=par(3); 
u=idct(uf); fu=lam+2*c2*u+3*c3*u.^2; p2pglob.fu=fu; 
Gu=speye(n); % dummy, returned here cause size used as dimension at places 