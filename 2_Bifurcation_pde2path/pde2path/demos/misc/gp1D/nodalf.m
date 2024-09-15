function f=nodalf(p,u) % nonlinearity, nodal version
par=u(p.nu+1:end); lam=par(1); s=par(2); del=par(3); sig=par(4); 
u=u(1:p.nu); x=getpte(p); x=x'; V=pot(x,s);
f=-V.*u-sig*u.^3-lam*u; 
