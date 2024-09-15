function r=slsG(p,u) % rhs for SL v_t=D*lap v-1/l-b*v+v^2/(1+v^2)
%                                 l_t=-D lap l+2cp*v+l*(rho+bp-2*v/(1+v^2)^2; 
par=u(p.nu+1:end); r=par(1); bp=par(2); cp=par(3); D=par(4); % extract pars 
v=u(1:p.np); l=u(p.np+1:2*p.np); % extract sol components 
f1=-1./l-bp*v+v.^2./(1+v.^2); % nonlin., first component 
f2=2*cp*v+l.*(r+bp-2*v./(1+v.^2).^2); % 2nd component 
f=[f1;f2]; 
r=D*p.mat.K*u(1:p.nu)-p.mat.M*f; % residual 