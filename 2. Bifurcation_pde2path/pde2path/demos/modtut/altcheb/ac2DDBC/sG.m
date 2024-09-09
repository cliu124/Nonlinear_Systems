function r=sG(p,u)  % PDE rhs 
n=p.nu; par=u(n+1:end); u=u(1:n); % split u into par and pde parts 
lam=par(2); c2=par(3); c3=par(4); f=lam*u+c2*u.^2+c3*u.^3; % nonlin. 
uf=zeros((p.nx+2)*(p.ny+2),1); uf(p.bui)=u(1:p.np); % full u 
r1=-par(1)*p.mat.L*uf;  r=r1(p.bui)-f;    % bulk part of PDE   