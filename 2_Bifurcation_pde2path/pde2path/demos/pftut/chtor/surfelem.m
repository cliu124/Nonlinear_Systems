function dS=surfelem(p) % for torus 
po=getpte(p); y=po(2,:)'; %y=(p.mat.drop*po(2,:)'); 
par=p.u(p.nu+1:end); R=par(4); rho=par(5); cy=cos(y); 
dS=(R+rho*cy)*rho; 
