function q=qf(p,u) % phase condition only for translation
q=p.u0x(1:p.nu)'*u(1:p.nu); 