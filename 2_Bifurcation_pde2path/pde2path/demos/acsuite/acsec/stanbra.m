function out=stanbra(p,u) % overload of stanbra to append the normalized L2-norm 
upde=u(1:p.nu); np=p.nu/p.nc.neq; 
l2=u(1:np)'*(p.mat.M(1:np,1:np)*u(1:np)); l2r=sqrt(l2/p.vol); 
out=[u(p.nu+1:end); max(abs(upde(1:np))); min(abs(upde(1:np)));l2r]; 