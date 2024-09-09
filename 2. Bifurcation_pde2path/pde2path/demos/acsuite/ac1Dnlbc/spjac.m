function Guuphi=spjac(p,u) % for ac1Dnlbc
par=u(2*p.nu+1:end); phi=u(p.nu+1:2*p.nu); u=u(1:p.nu); % params, Evec, PDE-vars
fuu=6*u-20*par(3)*u.^3; al=par(5); 
Guuphi=-p.mat.M*spdiags(fuu.*phi,0,p.nu,p.nu)...  % from the bulk
       +p.nc.sf*spdiags(2*al*p.mat.Q2*phi,0,p.nu,p.nu); % from the BC 