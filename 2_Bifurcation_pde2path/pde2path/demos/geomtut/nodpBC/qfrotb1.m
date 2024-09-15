function q=qfrot(p,u) 
% constraints: 3 translations, and rotation around z (for non S^1 solutions) 
u=u(1:p.nu); N0=getN(p,p.X); grX=grad(p.X,p.tri); % surface gradient
grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:); grXz=grX(2*p.nt+1:3*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; grXz=E*grXz; % interpol.to nodes 
dphix=grXx*p.X; dphix=dot(dphix,N0,2); 
dphiy=grXy*p.X; dphiy=dot(dphiy,N0,2);
dphiz=grXz*p.X; dphiz=dot(dphiz,N0,2);
qx=(p.mat.fill'*dphix)'*u; qy=(p.mat.fill'*dphiy)'*u; qz=(p.mat.fill'*dphiz)'*u;
phi=(-p.X(:,2).*grXx*p.X+p.X(:,1).*grXy*p.X)./(p.X(:,1).^2+p.X(:,2).^2); phi=dot(phi,N0,2);
qphi=(p.mat.fill'*phi)'*u; 



q=[qx;qy;qz;qphi];