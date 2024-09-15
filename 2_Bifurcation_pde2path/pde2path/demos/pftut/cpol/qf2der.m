function qu=qf2der(p,u) % derivative of smass-cons
par=u(p.nu+1:end); R=par(1); 
qmu=[R^2*p.mat.vM1, p.mat.vM2(1:end)]; 
qru=[(p.mat.Dphi*p.u(1:p.nus))', sparse(1,p.nu-p.nus)];  
qu=[qmu;qru]; 