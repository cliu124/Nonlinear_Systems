function q=qf(p,u)  % rotational PC 
Krot=p.mat.Dphi; uref=p.u(1:p.nu); 
q=(Krot*uref)'*(u(1:p.nu)-uref); 