function q=qV(p,u) % volume constraint 
par=u(p.nu+1:end); V0=par(2); u=u(1:p.nu); N0=getN(p,p.X); X=p.X+u.*N0;
N=cross(p.mat.Dx*X,p.mat.Dy*X,2); % normal at X, NOT normalized 
V=sum(p.mat.M*(dot(X,N,2)))/3; q=V-V0; 