function Gu=sGjac(p,u)  % AC Jacobian 
par=u(p.nu+1:end); lam=par(1); u=u(1:p.nu);  % split u into parameters and PDE variables 
K=p.mat.K;  fu=lam-3*u.^2; 
Fu=spdiags(fu,0,p.nu,p.nu);  % convert derivatives into (sparse) matrix 
Gu=K-p.mat.M*Fu+p.nc.sf*p.mat.Q;        % the Jacobian matrix 