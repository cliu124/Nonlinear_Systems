function qu=qjac(p,u) % q_u for rotational PC 
Krot=p.mat.Dphi; uref=p.u(1:p.nu); qu=(Krot*uref)'; 


