function Guuph=bpjac(p,u) % for BP and FP cont 
n=p.np; u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
ov=ones(n,1); den=(u2-u1).^2+1; d2=den.^2; d3=d2.*den; 
f1uu=-6*(u2-u1)./d2+8*(u2-u1).^3./d3; 
f1uv=-f1uu; f1vv=f1uu; 
f2uu=0*u1; f2uv=f1uv; f2vv=f1vv;
ph1=u(p.nu+1:p.nu+p.np); ph2=u(p.nu+p.np+1:2*p.nu);
M1=spdiags(f1uu.*ph1+f1uv.*ph2,0,n,n); 
M2=spdiags(f1uv.*ph1+f1vv.*ph2,0,n,n); 
M3=spdiags(f2uu.*ph1+f2uv.*ph2,0,n,n); 
M4=spdiags(f2uv.*ph1+f2vv.*ph2,0,n,n);
Guuph=-p.mat.M*[[M1 M2]; [M3 M4]]; 