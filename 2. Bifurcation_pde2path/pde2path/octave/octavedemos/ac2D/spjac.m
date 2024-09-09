function Guuph=spjac(p,u) 
% get partial_u (G_u phi), called in getGu if p.spcontsw==1 
par=u(2*p.nu+1:end); ph=u(p.nu+1:2*p.nu); u=u(1:p.nu); % params, Evec, PDE-vars
fuu=6*u-20*par(3)*u.^3; Guuph=-p.mat.M*spdiags(fuu.*ph,0,p.nu,p.nu); 
