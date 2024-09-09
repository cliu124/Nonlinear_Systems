function out=stanbra(p,u)
% mod of stanbra to also put the L2-norm on the branch 
uf=u(1:p.nu); np=p.nu/p.nc.neq; ux=p.mat.F'*uf; % inverse FFT 
ul2=sqrt(sum(ux.^2)/p.nu); 
out=[u(p.nu+1:end); max(abs(ux(1:np))); min(abs(ux(1:np))); ul2]; 