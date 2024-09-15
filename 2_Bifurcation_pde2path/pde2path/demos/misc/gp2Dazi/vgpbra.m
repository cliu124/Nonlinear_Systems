function out=vgpbra(p,u) % bd output for vGP
n=p.np; u1=u(1:n); u2=u(n+1:2*n); u3=u(2*n+1:3*n); u4=u(3*n+1:4*n); 
par=p.u(p.nu+1:end); ua1=u1.^2+u2.^2; ua2=u3.^2+u4.^2;
r1=max(abs(u1)); i1=max(abs(u2)); r2=max(abs(u3)); i2=max(abs(u4)); 
N1=sum(p.mat.M(1:n,1:n)*ua1,1);  N2=sum(p.mat.M(1:n,1:n)*ua2,1); N=N1+N2; 
out=[par; r1; i1; i1/r1; r2; i2; i2/r2; N1; N2; N]; 
%   1-5   6                             12  13  14