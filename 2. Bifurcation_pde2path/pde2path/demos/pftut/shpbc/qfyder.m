function qu=qfyder(p,u) % derivative of standard transl. phase condition
qu=(p.mat.Dy*p.u(1:p.nu))';  