function Guuph=bpjac(p,u) % for BP continuation in 'brussel' 
nu=p.nu; np=p.np; M=p.mat.M(1:np,1:np); 
u1=u(1:np); u2=u(np+1:2*np); 
ph1=u(nu+1:nu+np); ph2=u(nu+np+1:nu+2*np);  
Guuph=sparse(nu,nu); 
f1uu=2*u2; f1uv=2*u1; f2uu=-f1uu; f2uv=-f1uv; % other 2nd der.=0
M11=spdiags(f1uu.*ph1+f2uu.*ph2,0,np,np); 
M12=spdiags(f1uv.*ph1+f2uv.*ph2,0,np,np); 
M21=spdiags(f1uv.*ph1+f2uv.*ph2,0,np,np); 
Guuph(1:np,1:np)=-M*M11; Guuph(1:np,np+1:2*np)=-M*M12; 
Guuph(np+1:2*np,1:np)=-M*M21; 