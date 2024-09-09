function Gu=sGjac(p,u)  
par=u(p.nu+1:end); eps=par(2); u=u(1:p.nu);  % split u into parameters and PDE variables 
K=p.mat.K;  fu=1-3*u.^2; 
Fu=spdiags(fu,0,p.nu,p.nu);  % convert derivatives into (sparse) matrix 
Gu=eps^2*K-p.mat.M*Fu;        % the Jacobian matrix 