function r=sG(p,u)  % PDE rhs for AC
par=u(p.nu+1:end); lam=par(1); s=par(2); ga=par(3); u=u(1:p.nu); % split u into parameters and PDE variables 
f=lam*u-u.^3+ga; K=p.mat.K; Krot=p.mat.Krot; 
r=K*u-p.mat.M*f+s*Krot*u;    % the rhs 
