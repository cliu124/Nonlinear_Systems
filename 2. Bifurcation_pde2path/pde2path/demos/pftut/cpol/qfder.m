function qu=qfmder(p,u) % derivative of smass-cons
qu=[(p.mat.Dphi*p.u(1:p.nus))', sparse(1,p.nu-p.nus)];  %qu=p.u0x'; 