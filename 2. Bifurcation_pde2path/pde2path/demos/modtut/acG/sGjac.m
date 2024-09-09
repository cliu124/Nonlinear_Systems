function Gu=sGjac(p,u) % AC on graph, Jacobian 
n=p.nu; par=u(n+1:end);  u=u(1:n); lam=par(2); c2=par(3); c3=par(4); 
fu=lam+2*c2*u+3*c3*u.^2; Fu=spdiags(fu,0,n,n); % local Jac, turned into matrix 
Gu=par(1)*p.mat.L-Fu;  % Jac=c*Lap-f_u 