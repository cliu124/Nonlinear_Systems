function r=sG(p,u)  
par=u(p.nu+1:end); u=u(1:p.nu); f=par(2)*u+par(4)*u.^2-par(3)*u.^3; 
r=par(1)*p.mat.K*u-p.mat.M*f; 