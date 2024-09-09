function Guuphi=spjac(p,u) 
% partial_u (G_u phi), called in getGu if p.sw.spcont=1 and p.sw.spjac=1
par=u(2*p.nu+1:end); phi=u(p.nu+1:2*p.nu); u=u(1:p.nu); % params, Evec, PDE-vars 
fuu=6*u-20*par(3)*u.^3; Guuphi=-p.mat.M*spdiags(fuu.*phi,0,p.nu,p.nu); 
