function out=stanbra(p,u)
% mod of stanbra to also put the L2-norm on the branch 
u=u(1:p.nu); np=p.nu/p.nc.neq; 
ul2=sqrt(sum(u.^2)/p.nu); 
out=[u(p.nu+1:end); max(abs(u(1:np))); min(abs(u(1:np))); ul2]; 