function r=sG(p,u) % PDE rhs for CH; first split u into pars and PDE-vars 
par=u(p.nu+1:end); eps=par(2); lam=par(3); u=u(1:p.nu); 
f=u+lam-u.^3; r=eps^2*p.mat.K*u-p.mat.M*f;  %  nonlinearity and rhs 