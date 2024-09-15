function r=sG(p,u)  % ac2D 
par=u(p.nu+1:end); u=u(1:p.nu); R=par(1); lam=par(2); ga=par(3); s=par(4);  
f=lam*u+ga*u.^2-u.^3; 
r=p.mat.K*u/R^2-p.mat.M*f+s*p.mat.Dphi*u; 