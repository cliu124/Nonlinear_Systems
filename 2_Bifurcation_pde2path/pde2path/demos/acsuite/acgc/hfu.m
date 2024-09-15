function h=hfu(p,u) 
par=u(p.nu+1:end); j=par(5); u=u(1:p.nu); h=u.^j; 
