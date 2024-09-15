function r=sGenn(p,u)  % Enneper 
Xbc=bcX(p,u); p.X(p.idx,:)=Xbc; % set Bdcurve 
par=u(p.nu+1:end); H0=par(1); N=getN(p,p.X); X=p.X+u(1:p.np).*N; 
M=getM(p,X); LB=cotmatrix(X,p.tri); N=getN(p,X); 
r=-0.5*dot(LB*X,N,2)+M*(H0*ones(p.np,1)); % rhs-PDE, i.e., H-H0=0
r(p.idx)=u(p.idx); 