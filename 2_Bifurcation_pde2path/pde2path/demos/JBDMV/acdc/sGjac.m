function Gu=sGjac(p,u) % 'Jacobian' for AC dead core
a=afu(p,u); par=u(p.nu+1:end); lam=par(1); ga=par(2); c=par(3); n=p.np; 
u=u(1:p.nu); up=max(u,p.Jdel); % keeping u away from 0 for singular part 
fusing=ga*up.^(ga-1); % the singular part, 
fu=(lam-a)+2*(a+1).*u-3*u.^2-lam*fusing; Fu=spdiags(fu,0,n,n); 
Gu=c*p.mat.K-p.mat.M*Fu+p.nc.sf*p.mat.Q; % build Jac from bulk and BCs                