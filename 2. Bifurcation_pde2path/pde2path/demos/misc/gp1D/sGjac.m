function Gu=sGjac(p,u) % jac for real GP 
par=u(p.nu+1:end); lam=par(1); s=par(2); sig=par(3); ga=par(4); 
u=u(1:p.nu);  x=getpte(p); x=x'; V=pot(x,s);
fu=-V-ga*(2*sig+1)*u.^(2*sig)-lam; 
Gu=p.mat.K-p.mat.M*spdiags(fu,0,p.np,p.np); 