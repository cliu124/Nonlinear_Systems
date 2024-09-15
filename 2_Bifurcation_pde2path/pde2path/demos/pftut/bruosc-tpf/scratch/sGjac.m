function Gu=sGjac(p,u)
[f1u,f1v,f2u,f2v]=njac(p,u); n=(p.nu-2)/2; 
par=u(p.nu+1:end); lam=par(1); cp=par(3); del=par(4); om=par(5); 
u1=u(1:n); v1=u(p.nu-1); v2=u(p.nu);
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Ks=p.mat.K; K=[0*Ks -Ks; Ks p.mat.Ms]; 
Gu=K-p.mat.M0*Fu; 
Gv1=[cp*p.mat.Ms*u1; 0*u1]; 
av=[[del-3*v1^2-v2^2, -om-2*v1*v2]; [om-2*v1*v2, del-v1^2-3*v2^2]]; % osc.Jac 
Gu=[[Gu -Gv1 0*Gv1];[zeros(2,2*p.np) -av]]; % put it together
end 

function [f1u,f1v,f2u,f2v]=njac(p,u) 
n=(p.nu-2)/2; u1=u(1:n); v1=u(p.nu-1); 
par=u(p.nu+1:end); lam=par(1); nup=par(2); cp=par(3); ov=ones(n,1);  
if p.qc==1; f1u=-(1-lam-cp*v1)*ov+2*nup*u1-3*u1.^2; % quad-cub
else f1u=-(1-lam-cp*v1)*ov+3*nup*u1.^2-5*u1.^4; end 
f1v=-2*ov; f2u=0*ov; f2v=0*ov; 
end