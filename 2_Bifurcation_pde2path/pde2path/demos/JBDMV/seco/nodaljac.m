function [f1u,f1v,f2u,f2v]=nodaljac(p,u) % for seco 
n=p.np; u1=u(1:n); u2=u(n+1:2*n); den=(u2-u1).^2+1; % u, and a denominator
par=u(p.nu+1:end); al=par(2); tau=par(3); % parameters in nonlinearity 
f1u=-1./den+2*(u2-u1).^2./(den.^2)-tau; f1v=-f1u-tau; 
f2u=al*ones(n,1);  f2v=-f2u; 