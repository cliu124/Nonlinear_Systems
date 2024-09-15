function bc=hodummybc(ya,yb,p,par,T) 
% hodummybc: dummy BC for Hopf, used (and then irrelevant) in tomassempbc
bc=[ya(1:p.nu+p.nc.nq)]; 
