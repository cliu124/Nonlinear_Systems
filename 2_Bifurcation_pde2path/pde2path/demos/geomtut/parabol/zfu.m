function z=zfu(p,x,y) % set BCs here, use u=0 on bdry in sG
par=p.u(p.nu+1:end); a=par(1); b=par(2); 
z=x.^2/a^2+y.^2/b^2; 