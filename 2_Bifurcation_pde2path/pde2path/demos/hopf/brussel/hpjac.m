function Guuph=hpjac(p,u) % for HP continuation in 'brussel' 
nu=p.nu; np=p.np; M=p.mat.M(1:np,1:np); 
u1=u(1:np); u2=u(np+1:2*np); 
% u3=u(2*np+1:3*np); 
ph1r=u(nu+1:nu+np); ph2r=u(nu+np+1:nu+2*np);  
%ph3r=u(nu+2*np+1:nu+3*np);
ph1i=u(2*nu+1:2*nu+np); ph2i=u(2*nu+np+1:2*nu+2*np);  
%ph3i=u(2*nu+2*np+1:2*nu+3*np);

Guuph=sparse(2*nu,nu); 
f1uu=2*u2; f1uv=2*u1; 
f2uu=-f1uu; f2uv=-f1uv; 

M11=spdiags(f1uu.*ph1r+f1uv.*ph2r,0,np,np); M12=spdiags(f1uv.*ph1r,0,np,np); 
M21=spdiags(f2uu.*ph1r+f2uv.*ph2r,0,np,np); M22=spdiags(f2uv.*ph1r,0,np,np); 
Guuph(1:np,1:np)=-M*M11; Guuph(1:np,np+1:2*np)=-M*M12; 
Guuph(np+1:2*np,1:np)=-M*M21; Guuph(np+1:2*np,np+1:2*np)=-M*M22; 
M11=spdiags(f1uu.*ph1i+f1uv.*ph2i,0,np,np); M12=spdiags(f1uv.*ph1i,0,np,np); 
M21=spdiags(f2uu.*ph1i+f2uv.*ph2i,0,np,np); M22=spdiags(f2uv.*ph1i,0,np,np); 
Guuph(nu+1:nu+np,1:np)=-M*M11; Guuph(nu+1:nu+np,np+1:2*np)=-M*M12; 
Guuph(nu+np+1:nu+2*np,1:np)=-M*M21; Guuph(nu+np+1:nu+2*np,np+1:2*np)=-M*M22; 