function jc=polljcf(p,u) % current value for pollution 
par=u(p.nu+1:end); pr=par(2); vp=par(3); ga=par(5);
y=u(1:p.np); z=u(p.np+1:2*p.np);l1=u(2*p.np+1:3*p.np); % extract soln-components
k=-(1+l1)/ga; c=k+ga*k.^2/2; jc=pr*y-vp*z-c; % compute k, then J 