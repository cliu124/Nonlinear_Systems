function r=sG(p,u) % rhs for SH, see nodalf
f=nodalf(p,u); par=u(p.nu+1:end); Ks=p.mat.K; del=par(4); om=par(5); 
v1=u(p.nu-1); v2=u(p.nu); 
K=[0*Ks -Ks; Ks p.mat.Ms]; 
ru=K*u(1:p.nu-2)-p.mat.M0*f; 
va=v1^2+v2^2; rv1=-del*v1+om*v2+va*v1; rv2=-om*v1-del*v2+va*v2; % the oscill.eqns
r=[ru;rv1;rv2]; 