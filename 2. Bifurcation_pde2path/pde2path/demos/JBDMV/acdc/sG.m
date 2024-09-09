function r=sG(p,u)  % PDE rhs 
a=afu(p,u); par=u(p.nu+1:end); lam=par(1); ga=par(2); c=par(3); 
u=u(1:p.nu); up=max(u,0);   
fti=-u.*(u-1).*(u-a); f=fti-lam*(up.^ga-u); 
r=c*p.mat.K*u-p.mat.M*f+p.nc.sf*(p.mat.Q*u-p.mat.G); 