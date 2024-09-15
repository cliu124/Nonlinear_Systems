function qu=qfder(p,u) % derivative of standard transl. phase condition
qu=[(p.mat.Dx*p.u(1:p.nu))';(p.mat.Dy*p.u(1:p.nu))'];  %qu=p.u0x'; 