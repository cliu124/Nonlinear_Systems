function r=vegsG(p,u) % rhs for vegOC problem 
par=u(p.nu+1:end); rho=par(1); g=par(2); eta=par(3); % extract param 
d=par(4); del=par(5); beta=par(6); xi=par(7); rp=par(8); 
up=par(9); rw=par(10); pp=par(12); al=par(13); 
v=u(1:p.np); w=u(p.np+1:2*p.np); % extract soln-components, states 
l1=u(2*p.np+1:3*p.np); l2=u(3*p.np+1:4*p.np); % co-states lam1, lam2 
[e,h]=efu(p,u);  % get effort E and harvest h 
f1=(g*w.*v.^eta-d*(1+del*v)).*v-h; 
f2=rp*(beta+xi*v)-(up*v+rw).*w;  
f3=rho*l1-pp*al*h./v-l1.*(g*(eta+1)*w.*v.^eta-2*d*del*v-d-al*h./v)...
    -l2.*(rp*xi-up*w); 
f4=rho*l2-l1.*(g*v.^(eta+1))-l2.*(-up*v-rw); 
f=[f1;f2;f3;f4]; r=p.mat.K*u(1:p.nu)-p.mat.M*f; % the residual 