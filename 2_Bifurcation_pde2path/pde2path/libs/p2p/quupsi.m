function quups=quupsi(p,u,r) 
% quupsi: derivative of qu*psi, usually zero (for constraint q linear in u) 
quups=sparse(p.nc.nq,p.nu); 