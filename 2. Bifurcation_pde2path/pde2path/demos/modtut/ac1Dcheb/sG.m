function r=sG(p,u)  % PDE rhs 
n=p.nu; par=u(n+1:end); u=u(1:n); lam=par(2); c2=par(3); c3=par(4); % split u 
f=lam*u+c2*u.^2+c3*u.^3; f(1)=0; f(n)=0; % 'nonlinearity', zeroed on the bdry 
r=-par(1)*p.mat.L*u-f; % compute L*u on extended domain 
  