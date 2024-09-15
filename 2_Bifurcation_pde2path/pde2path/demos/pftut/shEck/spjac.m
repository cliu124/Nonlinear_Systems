function duGuph=spjac(p,u) % for quad-cub-SH 
n=p.np; u1=u(1:p.np); par=u(2*p.nu+1:end); nup=par(2); 
ov=ones(n,1); 
f1uu=2*nup*ov-6*u1; f1uv=0*ov; f1vv=f1uv; 
f2uu=f1uv; f2uv=f1uv; f2vv=f1uv;
ph1=u(p.nu+1:p.nu+p.np); ph2=u(p.nu+p.np+1:2*p.nu);
% general form, although only upper left block nonzero 
M1=spdiags(f1uu.*ph1+f1uv.*ph2,0,n,n); 
M2=spdiags(f1uv.*ph1+f1vv.*ph2,0,n,n); 
M3=spdiags(f2uu.*ph1+f2uv.*ph2,0,n,n); 
M4=spdiags(f2uv.*ph1+f2vv.*ph2,0,n,n);
duGuph=-p.mat.M*[[M1 M2]; [M3 M4]]; 