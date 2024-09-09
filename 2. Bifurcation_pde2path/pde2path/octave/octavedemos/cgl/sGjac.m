function Gu=sGjac(p,u)    % Jacobian for cGL
n=p.np; [f1u,f1v,f2u,f2v]=njac(p,u); % the specifics 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];  % put partial derivatives 
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]]; % into (sparse) F_u
Gu=p.mat.K-p.mat.M*Fu; % multiply by M and combine with Laplacian 
end 

function [f1u,f1v,f2u,f2v]=njac(p,u) % local (no spat.derivatives) Jacobian 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end);  % extract fields 
r=par(1); nu=par(2); mu=par(3); c3=par(4); c5=par(5);  % and parameters 
ua=u1.^2+u2.^2; 
f1u=r-2*u1.*(c3*u1-mu*u2)-c3*ua-4*c5*ua.*u1.^2-c5*ua.^2; 
f1v=-nu-2*u2.*(c3*u1-mu*u2)+mu*ua-4*c5*ua.*u1.*u2;
f2u=nu-2*u1.*(c3*u2+mu*u1)-mu*ua-4*c5*ua.*u1.*u2; 
f2v=r-2*u2.*(c3*u2+mu*u1)-c3*ua-4*c5*ua.*u2.^2-c5*ua.^2;  
end 