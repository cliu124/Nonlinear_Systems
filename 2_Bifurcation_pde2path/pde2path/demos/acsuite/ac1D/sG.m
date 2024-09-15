function r=sG(p,u)  % PDE rhs 
f=nodalf(p,u); % compute "nonlinearity" (everything but diffusion) 
par=u(p.nu+1:end); u=u(1:p.nu); % split in par and PDE u 
r=par(1)*p.mat.K*u-p.mat.M*f...   % bulk part of PDE 
  +p.nc.sf*(p.mat.Q*u-par(4)*p.mat.G); %  BC