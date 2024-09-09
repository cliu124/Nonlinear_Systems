function bc=hobc(ya,yb,p,par) 
% hobc: BC for Hopf, for tomsol via tom 
% y padded with aux vars for pBC 
nu=p.nu; phasebc=ya(1:nu)'*p.hopf.u0dot; 
bc=[ya(1:nu)-ya(nu+1:2*nu); yb(1:nu)-yb(nu+1:2*nu);  phasebc]; 
