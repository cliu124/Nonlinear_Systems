function Gu=gpsGjac(p,u)
n=p.np; par=u(p.nu+1:end); om=par(1); mu=par(2); pl=par(3); 
u1=u(1:n); u2=u(n+1:2*n); ua=u1.^2+u2.^2; dia=p.mat.pot-mu-ua; 
f1u=dia-2*u1.^2; f1v=-2*u1.*u2; 
f2u=-2*u1.*u2;    f2v=dia-2*u2.^2; 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Fu=-Fu; 
Kr=p.mat.Krot; Q=p.mat.Q; 
Gu=p.mat.K-p.mat.M*Fu+[0*Kr -om*Kr; om*Kr 0*Kr]+p.sf*[Q 0*Q; 0*Q Q]...
    +pl*[0*Q -speye(n); speye(n) 0*Q]; 

