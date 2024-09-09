function qu=qjacrot(p,u) % jac of 3 transl.constraints + 1 rotation (around z) 
grX=grad(p.X,p.tri); N0=getN(p,p.X);
grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:); grXz=grX(2*p.nt+1:3*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; grXz=E*grXz; 
dx=grXx*p.X; dx=dot(dx,N0,2); 
dy=grXy*p.X; dy=dot(dy,N0,2);
dz=grXz*p.X; dz=dot(dz,N0,2);
qux=(p.mat.fill'*dx)'; quy=(p.mat.fill'*dy)'; quz=(p.mat.fill'*dz)';
dphi=(-p.X(:,2).*grXx*p.X+p.X(:,1).*grXy*p.X)./(p.X(:,1).^2+p.X(:,2).^2); 
dphi=dot(dphi,N0,2); quphi=(p.mat.fill'*dphi)';
qu=[qux;quy;quz;quphi];