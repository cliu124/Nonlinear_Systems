function J=bpjac(p,u) % \pa_u (G_u' psi), called in getGu for BPcont p.sw.spcont=1 
psi=u(p.nu+1:2*p.nu); u=u(1:p.nu); 
fuu=-6*u; J=-spdiags(fuu.*psi,0,p.nu,p.nu)*p.mat.M';  
%J=-p.mat.M*spdiags(fuu.*psi,0,p.nu,p.nu); 