function r=sG(p,u)  % AC with x-dependent terms (K(x) from oosetfemops) 
par=u(p.nu+1:end);  u=u(1:p.nu); % split u into parameters and PDE vars 
x=getpte(p); x=x'; % extract the point coordinates from p and transpose 
f=par(2)*u+u.^3-par(3)*u.^5+0.5*x.*u; % f, with x-dependent term
r=p.mat.K*u-p.mat.M*f;  % bulk part of PDE 