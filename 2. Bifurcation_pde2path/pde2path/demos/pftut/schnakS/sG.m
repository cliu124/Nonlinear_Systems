function r=sG(p,u) % Schnakenberg on sphere
par=u(p.nu+1:end); R=par(4); Ks=p.mat.K/R^2; 
f=nodalf(p,u); s=par(5); K=[Ks 0*Ks; 0*Ks par(3)*Ks]; 
r=K*u(1:p.nu)-p.mat.M*f+s*p.mat.Dphi*u(1:p.nu); 