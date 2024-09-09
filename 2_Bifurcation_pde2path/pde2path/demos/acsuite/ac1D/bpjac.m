function J=bpjac(p,u) % \pa_u (G_u phi), called in getGubpco if p.sw.spjac=1 
par=u(2*p.nu+1:end); psi=u(p.nu+1:2*p.nu); u=u(1:p.nu); % pars, Evec, u  
fuu=6*u-20*par(3)*u.^3; J=-spdiags(fuu.*psi,0,p.nu,p.nu)*p.mat.M'; 