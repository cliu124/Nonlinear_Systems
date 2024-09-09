function r=sG(p,u)  % PDE rhs for CH 
par=u(p.nu+1:end); eps=par(2); lam=par(3); s=par(6); K=p.mat.K; 
u=u(1:p.nu); f=u+lam-u.^3; 
r=eps^2*K*u-p.mat.M*f+s*p.mat.Dphi*u(1:p.nu);  