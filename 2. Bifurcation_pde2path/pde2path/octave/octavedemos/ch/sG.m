function r=sG(p,u)  % PDE rhs for CH 
par=u(p.nu+1:end); eps=par(2); lam=par(3); 
u=u(1:p.nu); % split u into parameters and PDE variables 
f=u+lam-u.^3; K=p.mat.K; r=eps^2*K*u-p.mat.M*f;    % the rhs 
