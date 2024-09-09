function q=qfrot(p,u) 
% constraints: 3 translations, and rotation around z (for non S^1 solutions) 
u=u(1:p.nu); N0=getN(p,p.X); grX=grad(p.X,p.tri); % surface gradient
grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:); grXz=grX(2*p.nt+1:3*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; grXz=E*grXz; % interpol.to nodes 
dx=grXx*p.X; dx=dot(dx,N0,2); % x-translations
dy=grXy*p.X; dy=dot(dy,N0,2); % y-transl.
dz=grXz*p.X; dz=dot(dz,N0,2); % z-transl.
qx=(p.mat.fill'*dx)'*u; qy=(p.mat.fill'*dy)'*u; qz=(p.mat.fill'*dz)'*u;
phi=(-p.X(:,2).*grXx*p.X+p.X(:,1).*grXy*p.X)./(p.X(:,1).^2+p.X(:,2).^2); phi=dot(phi,N0,2);
qphi=(p.mat.fill'*phi)'*u; 
q=[qx;qy;qz;qphi];