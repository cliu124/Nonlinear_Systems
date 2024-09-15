function q=qfy(p,u) % standard phase condition for translational invariance in 'y' 
uold=p.u(1:p.nu); u0y=p.mat.Ky*uold; q=u0y'*u(1:p.nu);