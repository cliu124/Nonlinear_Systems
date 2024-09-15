function r=sG(p,u)  % PDE rhs 
n=p.nu; par=u(n+1:end); gbd=bdfu(p,u); u=u(1:n); % split u into par and pde parts 
lam=par(2); c2=par(3); c3=par(4); f=lam*u+c2*u.^2+c3*u.^3; % nonlin. 
uf=zeros(p.na*p.nr,1); uf(p.na+1:p.na+n)=u(1:n); % full u 
uf(1:p.na)=gbd; % fill in bd-values 
r1=-par(1)*p.mat.L*uf+par(5)*p.mat.Dphi*uf;  r=r1(p.na+1:p.na+n)-f;    % bulk part of PDE   