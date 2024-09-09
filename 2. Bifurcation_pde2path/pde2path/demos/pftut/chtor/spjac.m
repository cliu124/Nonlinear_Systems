function J=spjac(p,u) % \pa_u (G_u phi), getGu for foldcont (p.sw.spcont=2) 
phi=u(p.nu+1:2*p.nu); u=u(1:p.nu); 
fuu=-6*u; J=-p.mat.M*spdiags(fuu.*phi,0,p.nu,p.nu);  