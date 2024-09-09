function pot=pot(p) % potential 
x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)';
pot=x.^2+y.^2;
