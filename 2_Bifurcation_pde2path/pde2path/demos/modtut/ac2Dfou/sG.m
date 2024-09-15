function r=sG(p,u)  % PDE rhs 
n=p.nu; par=u(n+1:end); uf=u(1:n); lam=par(2); c2=par(3); c3=par(4); 
F=p.mat.F; u=F'*uf; f=lam*u+c2*u.^2+c3*u.^3; ff=F*f; % nonlin. 
r=par(1)*p.mat.L*uf-ff;    % resi 