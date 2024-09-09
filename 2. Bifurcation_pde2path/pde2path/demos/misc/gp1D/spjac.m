function Guuphi=spjac(p,u) % for 1D GP 
par=u(2*p.nu+1:end); sig=par(3); ga=par(4); 
phi=u(p.nu+1:2*p.nu); u=u(1:p.nu); 
fuu=-ga*(2*sig+1)*2*sig*u.^(2*sig-1); Guuphi=-p.mat.M*spdiags(fuu.*phi,0,p.nu,p.nu); 