function qu=qfder_sym(p,u) % phase condition jacs translation and rotation
qu=(p.mat.Kx*p.u(1:p.nu))'; qu=[qu; (p.mat.R*p.u(1:p.nu))'];