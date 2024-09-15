function bc=nlbc(p,u) % nonlin., x-dep. BC; s(x)=-(0.5+x+y), int_Ga s dGa <0 !
lam=u(p.nu+1); enum=max(p.mesh.e(5,:));
g=mat2str(0);q=[mat2str(lam) '*(0.5+x+y).*(1-u)']; 
bc=gnbcs(p.nc.neq,enum,q,g); 

