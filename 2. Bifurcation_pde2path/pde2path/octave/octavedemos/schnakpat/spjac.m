function Guuph=spjac(p,u)
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(2*p.nu+1:end); s=par(2); % sigma
n=p.np; ov=ones(n,1); 
f1uu=2*u2+2*s*ov; f1uv=2*u1+2*s*u2.^(-2); 
f1vv=-4*s*(u1-u2.^(-1)).*u2.^(-3)+2*s*u2.^(-4);
f2uu=-f1uu; f2uv=-f1uv; f2vv=-f1vv;
ph1=u(p.nu+1:p.nu+p.np); ph2=u(p.nu+p.np+1:2*p.nu);
M1=spdiags(f1uu.*ph1+f1uv.*ph2,0,n,n); 
M2=spdiags(f1uv.*ph1+f1vv.*ph2,0,n,n); 
M3=spdiags(f2uu.*ph1+f2uv.*ph2,0,n,n); 
M4=spdiags(f2uv.*ph1+f2vv.*ph2,0,n,n);
Guuph=-p.mat.M*[[M1 M2]; [M3 M4]]; 