function q=qx(p,u) % standard phase condition for translational invariance in 'x' 
u0x=p.mat.Dx(1:p.np,1:p.np)*p.u(1:p.np); q=u0x'*u(1:p.np);