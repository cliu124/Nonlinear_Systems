function qu=qfxder(p,u) % derivative of standard transl. phase condition
qu=(p.mat.Dx*p.u(1:p.nu))';  