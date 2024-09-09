function f=horhs(t,u,p,par) 
% horhs: rhs for Hopf, with aux vars for pBC in t
T=u(2*p.nu+1); u=[u(1:p.nu); par]; f=resi(p,u); fz=zeros(p.nu+1,1); 
f=[-T*f; fz]; 


