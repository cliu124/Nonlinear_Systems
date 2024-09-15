function r=sGhs(p,u)  % PDE rhs, discrete mean curvature 
par=u(p.nu+1:end); H0=par(1); sx=par(4); sy=par(5);
u=u(1:p.np); N=getN(p,p.X); X=p.X+u.*N; 
dphix=N(:,1); dphiy=N(:,2); % PC for cog (center of gravity)
M=getM(p,X); LB=cotmatrix(X,p.tri); N=getN(p,X); 
r=-0.5*dot(LB*X,N,2)+M*(H0*ones(p.np,1)+sx*dphix+sy*dphiy); % PDE 
grX=grad(X,p.tri); E=c2P(X,p.tri); 
grXz=grX(2*p.nt+1:3*p.nt,:); grXz=E*grXz; 
w=grXz*(X(:,1).^2+X(:,2).^2); %w=grXz*(X(:,1)); % also works 
r(p.idx)=w(p.idx); % Neumann BCs 