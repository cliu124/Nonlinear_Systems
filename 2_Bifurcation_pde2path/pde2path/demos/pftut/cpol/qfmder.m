function qu=qfmder(p,u) % derivative of smass-cons
par=u(p.nu+1:end); R=par(1); 
qu=[R^2*p.mat.vM1, p.mat.vM2(1:end)]; 
