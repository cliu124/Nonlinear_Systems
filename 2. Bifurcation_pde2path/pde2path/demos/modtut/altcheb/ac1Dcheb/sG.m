function r=sG(p,u)  % PDE rhs 
n=p.nu; par=u(n+1:end); u=u(1:n); % split u into par and pde parts 
lam=par(2); c2=par(3); c3=par(4); f=lam*u+c2*u.^2+c3*u.^3; % nonlin. 
u=[u(1);u;u(n)]; % add boundary points (according to NBCs) 
r1=-par(1)*p.mat.D2*u; % compute L*u on extended domain 
r=r1(2:n+1)-f;    % bulk part of PDE 
  