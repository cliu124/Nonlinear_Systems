function q=qf_sym(p,u) % phase condition for both translation and rotation
par=u(p.nu+1:end); uox=p.mat.Kx*p.u(1:p.nu); q=uox'*(u(1:p.nu)-p.u(1:p.nu)); 
q=[q; (p.mat.R*p.u(1:p.nu))'*(u(1:p.nu)-p.u(1:p.nu))+par(8)];