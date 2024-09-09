function quuph=spqjac(p,u) % \pa_u(q_u*phi), needed for fold continuation
quuph=sparse(p.nc.nq,p.nu);% here just 0-matrix 