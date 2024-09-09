function r=nlbsG(p,u)
lam=u(p.nu+1); kstar=u(p.nu+2:p.nu+3);ksq=kstar(1)^2+kstar(2)^2;
uf=p.mat.fill*u(1:p.nu);
ure=uf(1:p.np); uim=uf(p.np+1:2*p.np); ua=ure.^2+uim.^2;
f1=-p.sig*ua.*ure; f2=-p.sig*ua.*uim; 
f1=f1-(ksq-lam+p.mat.pot).*ure; f2=f2-(ksq-lam+p.mat.pot).*uim;
%build the nonlinear part of f and map it to the periodic domain
f=p.mat.drop*[f1;f2]; 
r=(p.mat.K-p.mat.Kadv)*u(1:p.nu)-p.mat.M*f; 
