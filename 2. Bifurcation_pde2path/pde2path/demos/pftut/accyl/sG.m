function r=sG(p,u) % rhs for 2 AC problems coupled via boundary 
np2=p.np2; par=u(p.nu+1:end); lam=par(2); c=par(1); R=par(4); s=par(5);  
u1=u(1:p.nu1); u2=u(p.nu1+1:p.nu1+np2); % the 2 fields 
f1=lam*u1+u1.^3-par(3)*u1.^5; % only 1 'nonlinearity' 
K=p.mat.Kphi/R^2+p.mat.Kz; % stiffness on cylinder 
r1=c*K*u1-p.mat.M1*f1+s*p.mat.Dphi*u1; 
r2=c*p.mat.K2*u2/R^2; % residuals 
% now set up coupling: sort u_1-u_2 on common bdry into residual r2
bcc=p.Q2*(p.S1*p.mat.fill*u1-p.S2*u2); 
r2=r2+p.sf*bcc; 
r=[r1;r2]; 