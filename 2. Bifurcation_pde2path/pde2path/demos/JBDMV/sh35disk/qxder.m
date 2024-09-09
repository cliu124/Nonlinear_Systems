function qu=qxder(p,u) % derivative of standard transl. phase condition in x
qx=(p.mat.Dx(1:p.np,1:p.np)*p.u(1:p.np))'; qu=[qx, 0*qx];   