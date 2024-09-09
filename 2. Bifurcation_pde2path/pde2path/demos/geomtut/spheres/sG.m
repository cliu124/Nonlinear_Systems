function r=sG(p,u)  % PDE rhs, discrete mean curvature 
n=p.np; par=u(p.nu+1:end); H=par(1); % K=par(3); sx=par(5); sy=par(6); sz=par(7); %V0=par(2);
u=u(1:n); N=getN(p,p.X); uN=u.*N; X=p.X+uN; M=getM(p); 
LB=cotmatrix(X,p.tri); N=getN(p,X); r=-0.5*dot(LB*X,N,2)+M*(H*ones(n,1)); 