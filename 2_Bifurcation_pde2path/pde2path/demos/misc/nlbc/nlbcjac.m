function bc=nlbcjac(p,u) % nonlin., x-dep. BC 
lam=u(p.nu+1); enum=max(p.mesh.e(5,:));
g=mat2str(0);qj=[mat2str(lam) '*(0.5+x+y).*(1-2*u)']; 
bc=gnbcs(p.nc.neq,enum,qj,g); 

