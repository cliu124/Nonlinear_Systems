function Gu=sGjac(p,u) 
par=u(p.nu+1:end); Ks=p.mat.K; 
K=[[par(1)*Ks par(2)*Ks];[par(3)*Ks Ks]]; 
n=p.np; [f1u,f1v,f2u,f2v]=njac(p,u); 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];  % put partial derivatives 
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]]; % into (sparse) Jac 
Gu=K-p.mat.M*Fu;                    % multiply by M and add Laplacian 
end 

function [f1u,f1v,f2u,f2v]=njac(p,u) % local (no spat.derivatives) Jacobian 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); % extract fields 
par=u(p.nu+1:end); a1=par(5); b1=par(6); % and pars 
f1u=a1-3*u1.^2+b1*u2;  f1v=b1*u1; 
f2u=-f1u; f2v=-f1v; 
end
