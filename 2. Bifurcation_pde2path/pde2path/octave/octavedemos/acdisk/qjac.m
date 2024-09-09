function qu=qjac(p,u) % q_u for rotational PC 
Krot=p.mat.Krot; uref=p.u(1:p.nu); qu=(Krot*uref)'; 


