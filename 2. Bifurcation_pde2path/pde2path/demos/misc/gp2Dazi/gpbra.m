function out=gpbra(p,u) % bd output for GP
np=p.np; u1=u(1:np); u2=u(np+1:2*np); par=p.u(p.nu+1:end); 
imax=max(u2); rmax=max(u1); ua=u1.^2+u2.^2;  N=sum(p.mat.M(1:np,1:np)*ua,1); 
out=[par; rmax; imax; imax/rmax; N]; 