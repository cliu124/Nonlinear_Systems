function q=qf(p,u) % standard phase condition for transl. invariance in 'x' 
uold=p.u(1:p.nu); u0x=p.mat.Kx*uold; q=u0x'*u(1:p.nu);