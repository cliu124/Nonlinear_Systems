function r=sG(p,u)  % ac2D,
f=nodalf(p,u);  % the nonlinearity
par=u(p.nu+1:end); u=u(1:p.nu); % split u into parameters and PDE-part 
r=par(1)*p.mat.K*u-p.mat.M*f... % the bulk part 
    +p.nc.sf*(p.mat.Q*u-par(4)*p.mat.G); % the boundary terms via stiff-spring