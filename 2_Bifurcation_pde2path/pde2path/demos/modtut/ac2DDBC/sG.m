function r=sG(p,u)  % PDE rhs 
n=p.nu; par=u(n+1:end); u=u(1:n); % split u into par and (active) field 
lam=par(2); c2=par(3); c3=par(4); f=lam*u+c2*u.^2+c3*u.^3; % nonlin. 
uf=zeros((p.nx+2)*(p.ny+2),1); % all u, to be filled by bulk and bdry values 
uf(p.bui)=u(1:p.np); % filling in bulk 
[xx,yy]=meshgrid(p.x,p.y); yy=yy(:); uf(p.rb)=par(5)*cos(pi/2*yy(p.rb)); % bd
r1=-par(1)*p.mat.L*uf;  % acting with L on full u, 
r=r1(p.bui)-f;          % extract active (bulk) DoFs 