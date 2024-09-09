function p=nlbpmm(p) % post mesh-modification procedure
p.mat.pot=pot(p); p.mat.poti=pdeintrp(p.mesh.p,p.mesh.t,p.mat.pot);
if p.sw.bcper>0; p=box2per(p); end 