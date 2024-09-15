function Guuph=hpjac(p,u) % for HP continuation in cGL
nu=p.nu; np=nu/2; M=p.mat.M(1:np,1:np); 
u1=u(1:np); u2=u(np+1:2*np); 
par=u(3*nu+1:end); mu=par(3); c3=par(4); c5=par(5); %par'
ua=u1.^2+u2.^2; 
ph1r=u(nu+1:nu+np); ph2r=u(nu+np+1:nu+2*np);  
ph1i=u(2*nu+1:2*nu+np); ph2i=u(2*nu+np+1:3*nu);  

f1uu=-2*(3*c3*u1-mu*u2)-12*c5*ua.*u1-8*c5*u1.^3; 
f1uv=2*mu*u1-2*c3*u2-8*c5*u2.*u1.^2-4*c5*ua.*u2; 
f1vv=-2*(c3*u1-mu*u2)+4*mu*u2-4*c5*ua.*u1-8*c5*u1.*u2.^2;  

f2uu=-2*(c3*u2+3*mu*u1)-4*c5*ua.*u2-8*c5*u1.^2.*u2; 
f2uv=-2*c3*u1-2*mu*u2-4*c5*ua.*u1-8*c5*u2.^2.*u1; 
f2vv=-2*(3*c3*u2+mu*u1)-12*c5*ua.*u2-8*c5*u2.^3; 

M11=M*spdiags(f1uu.*ph1r+f1uv.*ph2r,0,np,np); 
M12=M*spdiags(f1uv.*ph1r+f1vv.*ph2r,0,np,np); 
M21=M*spdiags(f2uu.*ph1r+f2uv.*ph2r,0,np,np); 
M22=M*spdiags(f2uv.*ph1r+f2vv.*ph2r,0,np,np); 
M31=M*spdiags(f1uu.*ph1i+f1uv.*ph2i,0,np,np); 
M32=M*spdiags(f1uv.*ph1i+f1vv.*ph2i,0,np,np); 
M41=M*spdiags(f2uu.*ph1i+f2uv.*ph2i,0,np,np); 
M42=M*spdiags(f2uv.*ph1i+f2vv.*ph2i,0,np,np); 
Guuph=-[M11 M12; M21 M22; M31 M32; M41 M42];  
% should be correct, and gives good cont, but hpjac doesn't fit 