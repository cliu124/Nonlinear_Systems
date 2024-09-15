function q=qf(p,u) % constraints: 3 translations
u=u(1:p.nu); N0=getN(p,p.X); grX=grad(p.X,p.tri); 
grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:); grXz=grX(2*p.nt+1:3*p.nt,:);
I=c2P(p.X,p.tri); % I interpolates from triangle centers to points
grXx=I*grXx; grXy=I*grXy; grXz=I*grXz; % grXx, grXy, grXz now act on nodal u  
dx=grXx*p.X; dx=dot(dx,N0,2); % x-translations
dy=grXy*p.X; dy=dot(dy,N0,2); % y-translations
dz=grXz*p.X; dz=dot(dz,N0,2); % z-translations
qx=(p.mat.fill'*dx)'*u; qy=(p.mat.fill'*dy)'*u; qz=(p.mat.fill'*dz)'*u;
q=[qx;qy;qz];