function out=stanbra(p,u)
% mod of stanbra, using p.nu1 instead of p.nu 
u1=p.mat.fill*u(1:p.nu1); l21=u1'*(p.mat.M0*u1); l21=sqrt(l21); 
out=[u(p.nu+1:end); l21; max(u(1:p.nu)); min(u(1:p.nu))]; 