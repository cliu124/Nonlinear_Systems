function V=getV(p,u)
par=u(p.nu+1:end); l=par(7); u=u(1:p.nu); N0=getN(p,p.X); 
X=p.X+u.*N0; M=massmatrix(X,p.tri,'full'); N=getN(p,X);
V=abs(sum(M*(dot(X,N,2)))/3);  V=V+pi*2*l/3; % circular BD addtion 