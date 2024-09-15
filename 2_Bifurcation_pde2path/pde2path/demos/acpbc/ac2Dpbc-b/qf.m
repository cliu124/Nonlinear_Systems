function q=qf(p,u) % standard phase condition for translational invariance in 'x' 
uold=p.u(1:p.nu); u=u(1:p.nu); u0x=p.mat.Kx*uold; u0y=p.mat.Ky*uold; 
q=[u0x'*u; u0y'*u]; 