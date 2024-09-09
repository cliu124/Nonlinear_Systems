function r=sG(p,u)  % generic PDE residual, sfem=+- 1 setting 
par=u(p.nu+1:end); u=u(1:p.nu); 
f=par(2)*u+u.^3-par(3)*u.^5; 
r=par(1)*p.mat.K*u-p.mat.M*f...   % bulk part of PDE 
  +p.nc.sf*(p.mat.Q*u-par(4)*p.mat.G); %  BC
end