function Gu=nlbsjac(p,u)
lam=u(p.nu+1);kstar=u(p.nu+2:p.nu+3); %load parameters
ksq=kstar(1)^2+kstar(2)^2;
uf=p.mat.fill*u(1:p.nu); 
ure=uf(1:p.np); uim=uf(p.np+1:2*p.np); ua=ure.^2+uim.^2;
uau=2*ure; uav=2*uim; 
f1u=-p.sig*(ua+uau.*ure)-(ksq+p.mat.pot-lam); 
f1v=-p.sig*uav.*ure; 
f2u=-p.sig*uau.*uim; 
f2v=-p.sig*(ua+uav.*uim)-(ksq+p.mat.pot-lam); 
n=p.np;
%the u-dependent part of the Jacobian  
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
%map to the periodic domain
Fud=p.mat.drop*Fu*p.mat.fill;
%build the full Jacobian
Gu=p.mat.K-p.mat.Kadv-p.mat.M*Fud; 
