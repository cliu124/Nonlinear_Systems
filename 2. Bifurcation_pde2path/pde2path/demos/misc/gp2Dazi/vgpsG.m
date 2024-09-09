function r=vgpsG(p,u) % 2-compo GP2D 
par=u(p.nu+1:end); u=u(1:p.nu); n=p.np; 
om=par(1); mu1=par(2); pl=par(3); mu2=par(5); % lagrange-pars for phase and rot  
if size(p.mat.pot,1)~=n; p.mat.pot=pot(p); end; 
u1=u(1:n); v1=u(n+1:2*n); u2=u(2*n+1:3*n); v2=u(3*n+1:4*n);
ua1=u1.^2+v1.^2; ua2=u2.^2+v2.^2; 
dia1=p.mat.pot-mu1-ua1; dia2=p.mat.pot-mu2-ua2; 
f1=(dia1+0.5*ua2).*u1; f2=(dia1+0.5*ua2).*v1; 
f3=(dia2+0.5*ua1).*u2; f4=(dia2+0.5*ua1).*v2; 
f=-[f1;f2;f3;f4]; 
K=p.mat.K; Kr=p.mat.Krot; Q=p.mat.Q; 
r=[K 0*K; 0*K K]*u-p.mat.M*f+om*[-Kr*v1;Kr*u1;-Kr*v2; Kr*u2] ...
    +p.sf*[Q*u1;Q*v1;Q*u2;Q*v2] ...
    + pl*[-v1; u1; -v2; u2];  % phase condition
