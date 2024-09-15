function Guuph=spjac(p,u) % for FP continuation in cGL
nu=p.nu; np=nu/2; u1=u(1:np); u2=u(np+1:2*np); 
par=u(2*nu+1:end); mu=par(3); c3=par(4); c5=par(5); ua=u1.^2+u2.^2; 
phi1=u(nu+1:nu+np); phi2=u(nu+np+1:nu+2*np);  

f1uu=-2*(3*c3*u1-mu*u2)-12*c5*ua.*u1-8*c5*u1.^3; 
f1uv=2*mu*u1-2*c3*u2-8*c5*u2.*u1.^2-4*c5*ua.*u2; 
f1vv=-2*(c3*u1-3*mu*u2)-4*c5*ua.*u1-8*c5*u1.*u2.^2; 

f2uu=-2*(c3*u2+3*mu*u1)-4*c5*ua.*u2-8*c5*u1.^2.*u2; 
f2uv=-2*c3*u1-2*mu*u2-4*c5*ua.*u1-8*c5*u2.^2.*u1; 
f2vv=-2*(3*c3*u2+mu*u1)-12*c5*ua.*u2-8*c5*u2.^3; 

M11=spdiags(f1uu.*phi1+f1uv.*phi2,0,np,np); 
M12=spdiags(f1uv.*phi1+f1vv.*phi2,0,np,np); 
M21=spdiags(f2uu.*phi1+f2uv.*phi2,0,np,np); 
M22=spdiags(f2uv.*phi1+f2vv.*phi2,0,np,np); 
Guuph=-p.mat.M*[M11 M12; M21 M22]; 