function q=qf(p,u) % phase condition only for translation
par=u(p.nu+1:end); q=p.u0x(1:p.nu)'*u(1:p.nu); %-p.u0(1:p.nu)); 