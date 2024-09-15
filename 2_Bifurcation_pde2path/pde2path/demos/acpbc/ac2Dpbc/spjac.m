function Guuphi=spjac(p,u) % pa_u (G_u phi), called in getGu if p.spcontsw==1 
nu=p.nu; n=p.np; par=u(2*nu+1:end); 
phip=u(nu+1:2*nu); up=u(1:nu); % params, Evec, PDE-vars (per.domain) 
phi=p.mat.fill*phip; u=p.mat.fill*up; % u and phi on extended domain
fuu=6*u-20*par(3)*u.^3; % 2nd derivative 
Guuphi=-p.mat.M0*spdiags(fuu.*phi,0,n,n)*p.mat.fill; % mapped back to per.domian

