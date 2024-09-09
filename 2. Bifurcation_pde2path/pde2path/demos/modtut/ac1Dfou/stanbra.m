function out=stanbra(p,u) % mod of stanbra to also put the L2-norm on the branch 
uf=u(1:p.nu); np=p.nu/p.nc.neq; upde=p.mat.F'*uf; ul2=sqrt(sum(upde.^2)/p.nu); 
out=[u(p.nu+1:end); max(abs(upde(1:np))); min(abs(upde(1:np))); ul2]; 