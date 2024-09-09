function Gu=sGjac(p,u) % jac for SH with GC 
global p2pglob; % p2pglob.avec, cvec computed here and used in, e.g., gclss
hj=hjac(p,u); p2pglob.avec=sum(p.mat.M*spdiags(hj,0,p.nu,p.nu))./p.Om; 
p2pglob.cvec=p.mat.M*fnl_a(p,u); 
[f1u,f1v,f2u,f2v]=njac(p,u); n=p.nu/2;
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=p.mat.K-p.mat.M*(Fu+fnljac(p,u)); 
end 

function [f1u,f1v,f2u,f2v]=njac(p,u) 
n=p.nu/2; u1=u(1:n); par=u(p.nu+1:end); lam=par(1); nup=par(2); ov=ones(n,1);  
f1u=-(1-lam)*ov+2*nup*u1-3*u1.^2; f1v=-2*ov;
f2u=0*ov; f2v=0*ov; 
end