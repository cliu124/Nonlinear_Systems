function Gu=sGjac(p,u)  % Jac for ac1Dnlbc (includes D_u of BC) 
par=u(p.nu+1:end); al=par(5); u=u(1:p.nu); 
fu=par(2)+3*u.^2-5*par(3)*u.^4; Fu=spdiags(fu,0,p.nu,p.nu);  
Gu=par(1)*p.mat.K-p.mat.M*Fu...
   +p.mat.Q1+p.nc.sf*spdiags(1+2*al*u,0,p.nu,p.nu)*p.mat.Q2; 