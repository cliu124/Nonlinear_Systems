function r=sG(p,u) % Schnakenberg on cone; 
par=u(p.nu+1:end); a=par(4); [Ks,~]=LBcone(p,a); eps=par(6); 
f=nodalf(p,u); K=eps^2*[Ks 0*Ks; 0*Ks par(3)*Ks]; 
r=K*u(1:p.nu)-p.mat.M*f; 