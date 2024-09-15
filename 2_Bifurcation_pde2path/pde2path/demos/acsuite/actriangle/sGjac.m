function Gu=sGjac(p,u)  % generic PDE Jacobian, sfem=+- 1 setting 
par=u(p.nu+1:end); u=u(1:p.nu); 
fu=par(2)+2*par(4)*u-3*par(3)*u.^2; 
Fu=spdiags(fu,0,p.nu,p.nu);  
Gu=par(1)*p.mat.K-p.mat.M*Fu; 