function q=qf(p,u) % standard phase condition for translational invariance in 'x' 
uold=p.u(1:p.nu); u0x=p.mat.Dphi*uold;  %u0x=p.u0x; 
q=u0x'*u(1:p.nu);