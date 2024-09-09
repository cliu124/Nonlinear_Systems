function r=sG(p,u) % AC with x-dependent terms; const.coeff. K in oosetfemops, 
% thus x-dependent diffusion implemented here 
par=u(p.nu+1:end); n=p.nu; u=u(1:n); % split u into parameters and PDE vars 
x=getpte(p); x=x'; % extract point coordinates from p 
f=par(2)*u+u.^3-par(3)*u.^5+0.5*x.*u; F=p.mat.M*f; % f, with x-dependent term, and F
r=(spdiags(1+0.1*x.^2,0,n,n)*p.mat.K-spdiags(0.2*x,0,p.np,p.np)*p.mat.Kx)*u-F;  