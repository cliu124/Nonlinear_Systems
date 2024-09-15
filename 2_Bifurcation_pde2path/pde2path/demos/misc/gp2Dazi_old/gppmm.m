function p=nlbpmm(p) % post mesh-modification procedure
p.mat.pot=pot(p); p.mat.poti=pdeintrp(p.mesh.p,p.mesh.t,p.mat.pot);
p=setfemops(p); 
end 
