function Gu=slsGjac(p,u) % Jacobian for SL 
par=u(p.nu+1:end); r=par(1); bp=par(2); cp=par(3); D=par(4); % extract pars  
v=u(1:p.np); q=u(p.np+1:2*p.np); % extract components 
aux=1+v.^2; hp=2*v./(aux.^2); hpp=(2-6*v.^2)./(aux.^3); 
f1p=-bp+hp; f1q=1./q.^2; f2p=2*cp-hpp.*q; f2q=r+bp-hp; n=p.np;
Fu=[[spdiags(f1p,0,n,n),spdiags(f1q,0,n,n)];
    [spdiags(f2p,0,n,n),spdiags(f2q,0,n,n)]];
Gu=D*p.mat.K-p.mat.M*Fu; 