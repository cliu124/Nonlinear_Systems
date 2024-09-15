function r=sG(p,u)  % ac2D 
par=u(p.nu+1:end); u=u(1:p.nu); f=par(2)*u+u.^3-par(3)*u.^5; 
r=par(1)*p.mat.K*u-p.mat.M*f; 