function q=qf(p,u) % standard phase condition for translational invariance in 'x' 
uold=p.u(1:p.nu1); uox=p.mat.Dphi*uold; q=uox'*u(1:p.nu1);