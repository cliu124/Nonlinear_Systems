function q=qf(p,u) % mass constraint int u1+u2 dx=0
M=p.mat.M(1:p.np,1:p.np); par=u(p.nu+1:end); m=par(4); 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); q=sum(M*(u1+u2))/p.vol-m;