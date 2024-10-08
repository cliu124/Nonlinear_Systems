function Gu=sGjac(p,u)
[f1u,f1v,f2u,f2v]=njac(p,u); n=p.nu/2;
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
K=p.mat.Ks; Ms=p.mat.M(1:p.np, 1:p.np); Ksys=[[0*K -K];[K Ms]];
Gu=Ksys-p.mat.M*Fu; 
end 

function [f1u,f1v,f2u,f2v]=njac(p,u) 
n=p.nu/2; u1=u(1:n);par=u(p.nu+1:end); lam=par(1); nup=par(2); ov=ones(n,1); 
f1u=(lam-1)*ov+2*nup*u1-3*u1.^2; f1v=-2*ov;
f2u=0*ov; f2v=0*ov; 
end