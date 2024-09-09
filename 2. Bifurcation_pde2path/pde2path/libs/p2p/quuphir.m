function quupr=quuphir(p,u,r) 
% quuphir: derivative of qu*phir, usually zero 
% (for constraint q linear in u) 
quupr=sparse(p.nc.nq,p.nu); 