function r=sG(p,u)  % real GP 
par=u(p.nu+1:end);  lam=par(1); s=par(2); sig=par(3); ga=par(4); 
u=u(1:p.nu); x=getpte(p); x=x'; V=pot(x,s);
f=-V.*u-ga*u.^(2*sig+1)-lam*u;  r=p.mat.K*u-p.mat.M*f; 