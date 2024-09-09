function Gu=sGjac(p,u)  % Jacobian for AC 
par=u(p.nu+1:end); u=u(1:p.nu);  % split u into parameters and PDE variables 
fu=par(2)+3*u.^2-5*par(3)*u.^4; % local derivative of 'nonlinearity'  
Fu=spdiags(fu,0,p.nu,p.nu);  % put derivatives into (sparse) matrix 
Gu=par(1)*p.mat.K-p.mat.M*Fu;    % the Jacobian matrix 