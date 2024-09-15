function q=qfx(p,u) % standard phase condition for translational invariance in 'x' 
uold=p.u(1:p.nu); u0x=p.mat.Kx*uold; q=u0x'*u(1:p.nu);