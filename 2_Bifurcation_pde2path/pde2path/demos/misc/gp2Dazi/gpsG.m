function r=gpsG(p,u) % GP2D 
par=u(p.nu+1:end); u=u(1:p.nu); np=p.np; 
om=par(1); mu=par(2); pl=par(3); rl=par(4); % lagrange-pars for phase and rot  
if size(p.mat.pot,1)~=np; p.mat.pot=pot(p); end; 
u1=u(1:np); u2=u(np+1:2*np);
ua=u1.^2+u2.^2; dia=p.mat.pot-mu-ua; 
f1=dia.*u1; f2=dia.*u2; f=-[f1;f2]; 
Kr=p.mat.Krot; Q=p.mat.Q; 
r=p.mat.K*u-p.mat.M*f+[-om*Kr*u2;om*Kr*u1]+p.sf*[Q*u1;Q*u2] ...
    +pl*[-u2; u1];  % phase condition
