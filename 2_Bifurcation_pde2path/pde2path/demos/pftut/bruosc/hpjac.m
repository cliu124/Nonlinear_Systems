function Guuph=hpjac(p,u) % for HP continuation in 'brussel' 
nu=p.nu; np=p.np; M=p.mat.M(1:np,1:np); 
u1=u(1:np); u2=u(np+1:2*np); 
ph1r=u(nu+1:nu+np); ph2r=u(nu+np+1:nu+2*np);  
ph1i=u(2*nu+1:2*nu+np); ph2i=u(2*nu+np+1:2*nu+2*np);  
f1uu=2*u2; f1uv=2*u1; f2uu=-f1uu; f2uv=-f1uv; 
f1vv=0*f1uu; f2vv=0*f1uu; 
M11=M*spdiags(f1uu.*ph1r+f1uv.*ph2r,0,np,np); 
M12=M*spdiags(f1uv.*ph1r+f1vv.*ph2r,0,np,np); 
M21=M*spdiags(f2uu.*ph1r+f2uv.*ph2r,0,np,np); 
M22=M*spdiags(f2uv.*ph1r+f2vv.*ph2r,0,np,np); 
M31=M*spdiags(f1uu.*ph1i+f1uv.*ph2i,0,np,np); 
M32=M*spdiags(f1uv.*ph1i+f1vv.*ph2i,0,np,np); 
M41=M*spdiags(f2uu.*ph1i+f2uv.*ph2i,0,np,np); 
M42=M*spdiags(f2uv.*ph1i+f2vv.*ph2i,0,np,np); 
Guuph=-[M11 M12; M21 M22; M31 M32; M41 M42];  
