function out=stanbra(p,u)
% mod of stanbra first doing inverse FFT 
uf=u(1:p.nu); np=p.nu/p.nc.neq; ux=p.mat.F'*uf; % inverse dct 
out=[u(p.nu+1:end); max(abs(ux(1:np))); min(abs(ux(1:np)))]; 