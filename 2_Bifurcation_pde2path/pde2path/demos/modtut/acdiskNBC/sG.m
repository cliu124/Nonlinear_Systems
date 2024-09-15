function r=sG(p,u)  % PDE rhs 
n=p.nu; par=u(n+1:end); u=u(1:n); % split u into par and pde parts 
lam=par(2); c2=par(3); c3=par(4); f=lam*u+c2*u.^2+c3*u.^3; % nonlin. 
f(1:p.na)=0; % zero out at bdry 
r=-par(1)*p.mat.L*u-f+par(5)*p.mat.Dphi*u;    % bulk part of PDE   