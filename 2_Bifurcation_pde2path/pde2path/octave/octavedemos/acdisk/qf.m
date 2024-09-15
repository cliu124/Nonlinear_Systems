function q=qf(p,u)  % rotational PC 
Krot=p.mat.Krot; uref=p.u(1:p.nu); 
q=(Krot*uref)'*(u(1:p.nu)-uref(1:p.nu)); 