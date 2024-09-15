function hj=hjac(p,u) 
par=u(p.nu+1:end); j=par(5); u=u(1:p.nu); hj=j*u.^(j-1); 
