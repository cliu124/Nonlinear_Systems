function Gu=sGjac(p,u) % jac for AC with x-depend. terms (K(x) from oosetfemops) 
par=u(p.nu+1:end); u=u(1:p.nu); 
x=getpte(p); x=x'; fu=par(2)+3*u.^2-5*par(3)*u.^4+0.5*x; 
Gu=p.mat.K-p.mat.M*spdiags(fu,0,p.np,p.np); 