function pot=pot(p) % periodic potential (as in the proceedings ENOC, 2014 )
x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)';
pot=exp(-x.^2).*cos(y);