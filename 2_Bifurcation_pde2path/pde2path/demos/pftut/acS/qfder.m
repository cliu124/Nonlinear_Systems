function qu=qfder(p,u) % derivative of standard transl. phase condition
qu=(p.mat.Dphi*p.u(1:p.nu))';  %qu=p.u0x'; 