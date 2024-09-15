function qu=qfder_trans(p,u) % phase condition jac for only translation
qu=(p.mat.Kx*p.u(1:p.nu))'; 
end

