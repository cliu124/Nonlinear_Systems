function q=fchqf(p,u)
par=u(p.nu+1:end); u=u(1:p.nu); 
q=p.eta*u/p.vol-par(4);
