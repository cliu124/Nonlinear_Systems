function qu=qfder(p,u) % derivative of standard transl. phase condition
qu=(p.mat.Kx*p.u(1:p.nu))';  