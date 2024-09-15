function Gu=sGjac(p,u)  % PDE Jacobian 
par=u(p.nu+1:end); eps=par(2); u=u(1:p.nu); % split u into pars and PDE vars
fu=1-3*u.^2; Fu=spdiags(fu,0,p.nu,p.nu); % f-deriv., put into sparse matrix 
Gu=eps^2*p.mat.K-p.mat.M*Fu;  % the Jacobian matrix 