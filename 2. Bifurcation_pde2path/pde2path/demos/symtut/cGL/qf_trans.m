function q=qf_trans(p,u) % phase condition only for translation
par=u(p.nu+1:end);   q=(p.mat.Kx*p.u(1:p.nu))'*(u(1:p.nu)-p.u(1:p.nu))+par(8); 