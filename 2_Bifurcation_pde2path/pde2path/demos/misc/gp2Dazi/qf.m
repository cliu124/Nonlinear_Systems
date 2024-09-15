function q=qf(p,u) % phase condition 
q=(p.mat.R*p.u(1:p.nu))'*(u(1:p.nu)-p.u(1:p.nu));