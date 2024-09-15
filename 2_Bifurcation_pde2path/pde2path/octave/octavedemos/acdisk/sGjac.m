function Gu=sGjac(p,u)  % AC Jacobian 
par=u(p.nu+1:end); lam=par(1); s=par(2); u=u(1:p.nu);  % split u into parameters and PDE variables 
K=p.mat.K;  Krot=p.mat.Krot; fu=lam-3*u.^2; 
Fu=spdiags(fu,0,p.nu,p.nu);  % convert derivatives into (sparse) matrix 
Gu=K+s*Krot-p.mat.M*Fu;        % the Jacobian matrix 