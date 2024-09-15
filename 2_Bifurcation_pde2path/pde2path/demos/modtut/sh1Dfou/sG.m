function r=sG(p,u)  % sh1D, rhs F-version 
n=p.nu; par=u(n+1:end); uf=u(1:n); lam=par(1); c2=par(2); c3=par(3); % split 
F=p.mat.F; u=F'*uf; ff=lam*u+c2*u.^2+c3*u.^3; f=F*ff; % "nonlinearity"  
r=p.mat.L*uf-f;     % residual 
  