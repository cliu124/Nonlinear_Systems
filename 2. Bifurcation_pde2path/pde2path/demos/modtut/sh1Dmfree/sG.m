function r=sG(p,u)  % matrix free, multipl.in p2pglob, needed in lssgmres anyway 
n=p.nu; par=u(n+1:end); uf=u(1:n);lam=par(1);c2=par(2);c3=par(3); 
u=idct(uf); ff=lam*u+c2*u.^2+c3*u.^3; f=dct(ff); % "nonlinearity" via dct (no F)  
global p2pglob; r=p2pglob.mu.*uf-f;  % residual 
  