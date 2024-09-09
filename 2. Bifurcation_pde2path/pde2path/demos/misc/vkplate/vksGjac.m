function Gu=vksGjac(p,u)  % sfem=1 jacobian for von karman-plate
lam=u(p.nu+1); n=p.np; u5=u(4*n+1:5*n); u6=u(5*n+1:6*n); 
u7=u(6*n+1:7*n); u8=u(7*n+1:8*n); u9=u(8*n+1:9*n); u10=u(9*n+1:10*n); 
f2u5=spdiags(-u9,0,n,n); f2u6=spdiags(-u8,0,n,n); 
f2u7=spdiags(2*u10,0,n,n); f2u8=spdiags(-u6,0,n,n); 
f2u9=spdiags(-u5,0,n,n); f2u10=spdiags(2*u7,0,n,n);
f4u5=spdiags(u6,0,n,n); f4u6=spdiags(u5,0,n,n); 
f4u7=spdiags(-2*u7,0,n,n); zd=spdiags(zeros(n,1),0,n,n); 
Fu=[sparse([],[],[],n,10*n,0); 
    [zd zd zd zd f2u5 f2u6 f2u7  f2u8  f2u9 f2u10]; 
    sparse([],[],[],n,10*n,0); 
    [zd zd zd zd f4u5 f4u6 f4u7   zd    zd   zd]; 
    sparse([],[],[],6*n,10*n,0)]; 
Gu=p.mat.K+lam*p.mat.K2-p.mat.M*Fu;
end 
 