function q=qfy(p,u) % standard phase condition for translational invariance in 'x' 
uold=p.u(1:p.nu); u0y=p.mat.Dy*uold; q=u0y'*u(1:p.nu);