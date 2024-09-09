function r=vsG(p,u)  % complex GP 
par=u(p.nu+1:end);  lam=par(1); s=par(2); sig=par(3); ga=par(4); np=p.np; 
u1=u(1:np); u2=u(np+1:2*np); ua=u1.^2+u2.^2; 
x=getpte(p); x=x'; V=pot(x,s);
f1=-V.*u2-ga*ua.^sig.*u2-lam*u2; 
f2=V.*u1+ga*ua.^sig.*u1+lam*u1; 
f=[f1;f2];  r=p.mat.K*[u1;u2]-p.mat.M*f; 