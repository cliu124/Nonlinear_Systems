function Gu=sGjac(p,u) % Jacobian for cGL 
par=u(p.nu+1:end); del=par(7); n=p.np; [f1u,f1v,f2u,f2v]=njac(p,u); 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];  % put partial derivatives 
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]]*p.mat.fill; % into (sparse) Jac 
Gu=del^2*p.mat.K-p.mat.M0*Fu-par(6)*p.mat.Kx; % multiply by M and add Laplacian 