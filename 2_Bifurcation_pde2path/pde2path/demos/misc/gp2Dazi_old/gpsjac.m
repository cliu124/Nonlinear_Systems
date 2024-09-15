function Gu=gpsjac(p,u)
mu=u(p.nu+1); ga=u(p.nu+2); u=p.mat.fill*u(1:p.nu); 
fu=-p.mat.pot+mu-3*ga*u.^2; 
fu=spdiags(fu,0,p.nu,p.nu);
Gu=p.mat.K-p.mat.M*fu; 

