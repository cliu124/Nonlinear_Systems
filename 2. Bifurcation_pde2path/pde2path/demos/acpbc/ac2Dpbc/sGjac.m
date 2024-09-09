function Gu=sGjac(p,u)  % PDE Jacobian for AC with pBC 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
u=p.mat.fill*up; % extend ('fill') u to full domain 
fu=par(2)+3*u.^2-5*par(3)*u.^4; % local derivative on extended domain 
Fu=p.mat.M0*(spdiags(fu,0,p.np,p.np)*p.mat.fill); % map fu to per.dom
Gu=par(1)*p.mat.K+p.nc.sf*p.mat.Q-Fu; 