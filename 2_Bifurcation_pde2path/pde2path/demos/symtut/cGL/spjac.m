function Gvvph=spjac(p,u)  % for fold-cont in cGL, 
nu=p.nu; n=nu/2; par=u(2*nu+1:end); 
%u=p.mat.fill*u(1:p.nu); simplified version, omit fill for spjac 
u1=u(1:n); u2=u(n+1:2*n); mu=par(4); c3=par(5); c5=par(6); % sol and param.
ua=u1.^2+u2.^2; 
f1uu=-2*(c3*u1-mu*u2)-4*c3*u1-8*c5*u1.*(ua+u1.^2)-4*c5*ua.*u1; 
f1uv=2*mu*u1-2*c3*u2-4*c5*u2.*(2*u1.^2+ua); 
f1vv=-2*(c3*u1-mu*u2)+4*mu*u2-4*c5*(ua.*u1+2*u2.^2.*u1); 
f2uu=-2*(c3*u2+mu*u1)-4*mu*u1-4*c5*(ua.*u2+2*u1.^2.*u2); 
f2uv=-2*c3*u1-2*mu*u2-4*c5*(ua.*u1+2*u2.^2.*u1); 
f2vv=-2*(c3*u2+mu*u1)-4*c3*u2-8*c5*u2.*(ua+u2.^2)-4*c5*ua.*u2; 
ph1=u(2*n+1:3*n); ph2=u(3*n+1:4*n);
M1=spdiags(f1uu.*ph1+f1uv.*ph2,0,n,n); 
M2=spdiags(f1uv.*ph1+f1vv.*ph2,0,n,n); 
M3=spdiags(f2uu.*ph1+f2uv.*ph2,0,n,n); 
M4=spdiags(f2uv.*ph1+f2vv.*ph2,0,n,n);
Gvvph=-p.mat.M*[[M1 M2]; [M3 M4]]; 