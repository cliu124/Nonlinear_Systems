function Gu=sGjac(p,u)  % AC Jacobian 
par=u(p.nu+1:end); u=u(1:p.nu);  % split u into parameters and PDE variables 
K=par(1)*p.mat.K; fu=par(2)+3*u.^2-5*par(3)*u.^4; % K, and local derivatives 
Fu=spdiags(fu,0,p.nu,p.nu);  % put derivatives into (sparse) matrix 
Gu=K-p.mat.M*Fu;        % the Jacobian matrix 