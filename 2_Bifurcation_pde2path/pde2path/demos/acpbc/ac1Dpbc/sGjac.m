function Gu=sGjac(p,u)  % PDE Jacobian for AC with pBC 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
u=p.mat.fill*up; % extend ('fill') u to full domain 
x=getpte(p); x=x'; fu=par(2)+3*u.^2-5*par(3)*u.^4+0.5*x; % f_u on ext. domain 
Fu=p.mat.M0*(spdiags(fu,0,p.np,p.np)*p.mat.fill); % map fu to per.dom
Gu=p.mat.K-Fu; 