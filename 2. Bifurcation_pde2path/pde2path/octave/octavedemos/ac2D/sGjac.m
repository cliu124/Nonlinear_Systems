function Gu=sGjac(p,u)  % generic PDE Jacobian, sfem=+- 1 setting 
par=u(p.nu+1:end); u=u(1:p.nu); 
fu=par(2)+3*u.^2-5*par(3)*u.^4; 
Fu=spdiags(fu,0,p.nu,p.nu);  
Gu=par(1)*p.mat.K+p.nc.sf*p.mat.Q-p.mat.M*Fu; 