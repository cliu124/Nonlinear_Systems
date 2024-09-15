function Gu=sGjac(p,u)  % generic PDE Jacobian, sfem=+- 1 setting 
par=u(p.nu+1:end); up=u(1:p.nu); u=p.mat.fill*up; 
po=getpte(p);  r2=(po(1,:)+1).^2+(po(2,:)+1).^2+(po(3,:)+1).^2; r2=r2'; 
fu=par(2)+3*u.^2-5*par(3)*u.^4+1./(1+r2); 
Fu=p.mat.M0*(spdiags(fu,0,p.np,p.np)*p.mat.fill);  
Gu=par(1)*p.mat.K-Fu; %+p.nc.sf*p.mat.Q; 