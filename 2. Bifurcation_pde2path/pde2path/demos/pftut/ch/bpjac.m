function J=bpjac(p,u) % \pa_u (G_u' psi) for BPcont 
psi=u(p.nu+1:2*p.nu); u=u(1:p.nu); % adjoint kernel-vector and PDE-u
fuu=-6*u; J=-spdiags(fuu.*psi,0,p.nu,p.nu)*p.mat.M';  