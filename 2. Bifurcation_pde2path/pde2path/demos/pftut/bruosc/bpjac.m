function Guuph=bpjac(p,u)
nu=p.nu; np=p.np; M=p.mat.M(1:np,1:np); 
u1=u(1:np); u2=u(np+1:2*np); 
ph1=M*u(nu+1:nu+np); ph2=M*u(nu+np+1:nu+2*np);  
f1uu=2*u2; f1uv=2*u1; f2uu=-f1uu; f2uv=-f1uv; 
M11=spdiags(f1uu.*ph1+f2uu.*ph2,0,np,np); 
M12=spdiags(f1uv.*ph1+f2uv.*ph2,0,np,np); 
M21=spdiags(f1uv.*ph1+f2uv.*ph2,0,np,np); 
%Guuph=sparse(nu,nu); Guuph(1:np,1:np)=-M11; Guuph(1:np,np+1:2*np)=-M12; Guuph(np+1:2*np,1:np)=-M21;
Guuph=-[[M11, M12];[M21, 0*M21]];  