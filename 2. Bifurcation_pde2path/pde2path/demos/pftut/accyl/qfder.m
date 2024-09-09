function qu=qfder(p,u) % derivative of standard transl. phase condition
qu=[(p.mat.Dphi*p.u(1:p.nu1))' zeros(1,p.np2)]; 