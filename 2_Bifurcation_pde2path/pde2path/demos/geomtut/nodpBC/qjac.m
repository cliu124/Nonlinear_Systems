function qu=qjac(p,u) % jac of 3 transl.constraints 
grX=grad(p.X,p.tri); % surface grad; now extract components and interpol. to nodes
grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:); grXz=grX(2*p.nt+1:3*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; grXz=E*grXz; N0=getN(p,p.X);
dx=grXx*p.X; dx=dot(dx,N0,2); 
dy=grXy*p.X; dy=dot(dy,N0,2);
dz=grXz*p.X; dz=dot(dz,N0,2);
qux=(p.mat.fill'*dx)'; quy=(p.mat.fill'*dy)'; quz=(p.mat.fill'*dz)';
qu=[qux;quy;quz];