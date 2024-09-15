function f=nodalf(p,u) % for pollution 
par=u(p.nu+1:end); rho=par(1); pr=par(2); beta=par(3); ap=par(4); ga=par(5);
y=u(1:p.np); z=u(p.np+1:2*p.np); % extract soln-components, states 
l1=u(2*p.np+1:3*p.np); l2=u(3*p.np+1:4*p.np); % co-states 
k=-(1+l1)./ga; % control 
f1=-k; f2=y-z.*(1-z); f3=rho*l1-pr-l2; f4=(rho+1-2*z).*l2+beta; 
f=[f1; f2; f3; f4]; % 'nonlinearity' 