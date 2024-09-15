function r=sG(p,u)  % AC on graph rhs 
par=u(p.nu+1:end); u=u(1:p.nu); % split in par and PDE u 
lam=par(2); c2=par(3); c3=par(4); f=lam*u+c2*u.^2+c3*u.^3; % "nonlinearity"
r=par(1)*p.mat.L*u-f; % residual 