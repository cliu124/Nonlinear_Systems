function r=sGsc(p,u)  % spherical cap PDE (more generally: any CMC with DBCs) 
par=u(p.nu+1:end); H0=par(1); u=u(1:p.np); % split into PDE-u and parameters
N0=getN(p,p.X); X=p.X+u.*N0; N=getN(p,X); % base normal, new X, new normal   
M=getM(p,X); LB=cotmatrix(X,p.tri);  % mass matrix and Laplace-Beltrami
r=-0.5*dot(LB*X,N,2)+M*(H0*ones(p.np,1)); % rhs-PDE, i.e., -H(X)+H0=0
r(p.idx)=u(p.idx);    % Dirichlet BCs 