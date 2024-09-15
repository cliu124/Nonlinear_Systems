function V=getV(p,u)
par=u(p.nu+1:end); l=par(6); up=u(1:p.nu); u=p.mat.fill*up; N0=getN(p,p.X); 
X=p.X+u.*N0; M=getMf(p,X); N=getN(p,X);
V=2*abs(sum(M*(dot(X,N,2)))/3);  V=V+pi*2*l/3; % add bottom and top, why l? 