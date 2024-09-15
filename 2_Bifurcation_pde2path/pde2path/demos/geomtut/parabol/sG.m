function r=sG(p,u)  % PDE rhs 
par=u(p.nu+1:end); a=par(1); b=par(2);  
N0=getN(p,p.X); u1=u(1:p.nu); X=p.X+u1.*N0; % 
N=getN(p,X); M=getM(p,X); LB=cotmatrix(X,p.tri); 
x=X(:,1); y=X(:,2); a2=a^2; b2=b^2; a4=a2^2; b4=b2^2; 
Hf=(a2+b2+4*x.^2/a2+4*y.^2/b2)./(a2*b2*(1+4*x.^2/a4+4*y.^2/b4).^1.5); 
r=-0.5*dot(LB*X,N,2)+M*Hf; 
r(p.idx)=u(p.idx); 
