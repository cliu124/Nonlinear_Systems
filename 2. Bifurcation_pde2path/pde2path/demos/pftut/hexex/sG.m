function r=sG(p,u)  % generic PDE residual, sfem=+- 1 setting 
par=u(p.nu+1:end); u=u(1:p.nu); 
f=par(1)*(u+u.^3); 
r=p.mat.K*u-p.mat.M*f...   % bulk part of PDE 
  +p.nc.sf*p.mat.Q*u; %  BC
end