function r=sG(p,u) % Schnakenberg on torus 
Ks=p.mat.K; par=u(p.nu+1:end); % lam,sig,d,R,rho,s, now check R or rho are active
if any(ismember(p.nc.ilam,[4 5])); R=par(4); rho=par(5); Ks=LBtor(p,R,rho); end  
f=nodalf(p,u); s=par(6); K=[Ks 0*Ks; 0*Ks par(3)*Ks]; 
r=K*u(1:p.nu)-p.mat.M*f+s*p.mat.Dphi*u(1:p.nu); 